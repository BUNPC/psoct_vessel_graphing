function Mosaic3D_Telesto(ParameterFile, modality, nWorker)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIAL SETTING 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/autofs/space/vault_021/users/dv896/git_repo/projects/chi_real_data_code'));
addpath(genpath('/autofs/space/vault_021/users/dv896/git_repo/projects/brainreg'));
addpath(genpath('/autofs/space/vault_021/users/dv896/git_repo/projects/chi_real_data_code'));
% poolobj = parpool(10);

load(ParameterFile); %Parameters Mosaic3D Scan


%%% display numbers of workers 
fprintf('--nWorker is %i--\n', nWorker);

%%% display modality 
fprintf('--Modality is %s--\n', modality);

%%%
UseGPUflag = true;
if UseGPUflag; fprintf('--GPU option is considered--\n');
else; fprintf('--GPU option is not considered--\n');end

if isdeployed
    disp('--App is deployed--');
else
    disp('--App is NOT deployed--');
    addpath(genpath('/autofs/cluster/octdata2/users/Hui/tools/rob_utils/OCTBasic'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Begin!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING SLICE INPUT 
sliceidx    = Mosaic3D.sliceidx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING MOSAIC PARAMETERS 

fprintf(' - Loading Experiment file...\n %s\n', Mosaic3D.Exp);
S = whos('-file',Mosaic3D.Exp);
if any(contains({S(:).name},'Experiment_Fiji'))
    idx = find(contains({S(:).name},'Experiment_Fiji'));    isFiji = true;
elseif any(contains({S(:).name},'Experiment'))
    idx = find(contains({S(:).name},'Experiment'));         isFiji = false;
end
load(Mosaic3D.Exp,S(idx).name);
Experiment  = eval(S(idx).name);
fprintf(' - %s is loaded ...\n %s\n', S(idx).name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING DIRECTORIES 
indir       = Mosaic3D.indir;
outdir      = Mosaic3D.outdir; if ~exist(outdir,'dir'); mkdir(outdir); end
filetype    = Mosaic3D.InFileType; % 'nifti';

outdir = '/autofs/cluster/octdata2/users/Hui/20220201_SP002/process_run1//mus/mm_despeckle_dBI180';
mkdir(outdir)

fprintf(' - Input directory = \n%s\n',indir{:});
fprintf(' - Output directory = \n%s\n',outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ADDING CLIPPING
switch Scan.System
    case 'Octopus'
        XPixClip    = Parameters.YPixClip; %18;
        YPixClip    = Parameters.XPixClip; %0;
    case 'Telesto'
        XPixClip    = Parameters.XPixClip; %18;
        YPixClip    = Parameters.YPixClip; %0;
end
fprintf('XPixClip = %i\nYPixClip = %i\nSystem is %s\n\n', XPixClip, YPixClip,Scan.System);



% Set mosaic params
MapIndex = Experiment.MapIndex_Tot_offset+ Experiment.First_Tile - 1;

% Set unused tiles to -1 in MapIndex
% MapIndex(Experiment_Fiji.MapIndex_Fiji==-1)=-1;

% Set other parameters
MZL     = 180; % Mosaic3D.MZL;
fprintf('save %ipx of orignal %ipx\n', MZL, Scan.CropDepth);

sxyz    = Mosaic3D.sxyz;
%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 
NbPix = Experiment.NbPix;

X       = Experiment.X_Mean;           Y       = Experiment.Y_Mean;         % Y = fliplr(Y);
X       = X-min(X(:))+1;               Y       = Y-min(Y(:))+1;
sizerow = NbPix-XPixClip;              sizecol = NbPix-YPixClip;
MXL     = max(X(:))+sizerow-1;         MYL     = max(Y(:))+sizecol-1;
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% GENERATE RAMP

ramp_x = sizerow-round(median(diff(Experiment.X_Mean,[],1),'all','omitnan'));
ramp_y = sizecol-round(median(diff(Experiment.Y_Mean,[],2),'all','omitnan'));

ramp_xv       = 0:ramp_x-1;            ramp_yv       = 0:ramp_y-1;
x             = ones(1,sizerow);       y             = ones(1,sizecol);
x(1+ramp_xv)  = mat2gray(ramp_xv);     y(1+ramp_yv)  = mat2gray(ramp_yv);
x(end-ramp_xv)= mat2gray(ramp_xv);     y(end-ramp_yv)= mat2gray(ramp_yv);

RampOrig=x.'*y;
RampOrig3 = repmat(RampOrig, [1,1,MZL]);
% %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


tabrow=size(Experiment.MapIndex_Tot,1);
tabcol=size(Experiment.MapIndex_Tot,2);

% Loop slices
%poolobj = parpool(nWorker);
%parfor s = 1:size(sliceidx,2)
for s = 1:size(sliceidx,2)
    tic0=tic;
            
    sliceid_in  = sliceidx(1,s); 
    sliceid_out = sliceidx(2,s);
    sliceid_run = sliceidx(3,s); indir_curr = indir{sliceid_run};
    
%     if isFiji
%         MapIndex = squeeze(Experiment.MapIndex_Tot(:,:,sliceid_out));
%     else
    MapIndex = Experiment.MapIndex_Tot_offset+Experiment.First_Tile -1;
%     end
    % Start mosaic
    fprintf('\nStarting %s mosaic slice %d from run %d slice %d ...\n',modality,sliceid_out, sliceid_run, sliceid_in);

    
    M  = zeros(MXL,MYL,MZL,'single'); %flipped xy for telesto
    Ma = zeros(MXL,MYL,'single');
    Mosaic = zeros(MXL,MYL,MZL,'single');
    Masque = zeros(MXL,MYL,'single');
    
    for ii=1:tabcol
        fprintf('col %d / %d\n',ii,tabcol)
        for jj=1:tabrow
%             fprintf('row %d / %d\n',jj,tabrow)

            if MapIndex(jj,ii)>0 && ~isnan(X(jj,ii))
                columns =Y(jj,ii):Y(jj,ii)+sizecol-1;
                row     =X(jj,ii):X(jj,ii)+sizerow-1;
                
                currtile = (sliceid_in-1)*Experiment.TilesPerSlice+MapIndex(jj,ii);
                fprintf('currtile = %i\n',currtile);
                switch filetype
%                     case 'mat'
%                         Imag = [indir_curr filesep modality '_' sprintf('%03i',currtile) '.mat'];
%                         S = whos('-file',Imag);
%                         load(Imag);
%                         I = eval(S.name);
%                         if canUseGPU(); I =  gpuArray(I);end
%                         data = I(XPixClip+1:end,YPixClip+1:end,1:MZL);

                    case 'mat'
                        ftile = [indir_curr filesep modality '_' sprintf('%03i',currtile) '.mat'];
                        S = whos('-file',ftile);
                        load(ftile);
                        Imag = eval(S.name);
                    case 'nifti'
                        ftile = [indir_curr filesep 'test_processed_' sprintf('%03i',currtile) '_cropped.nii'];
                        Imag =  readnifti(ftile);
                end  
                
                Imag = single(Imag);

                if [canUseGPU() & UseGPUflag]; Imag =  gpuArray(Imag);end

                if strcmpi(Scan.System, 'Octopus');Imag = permute(Imag,[2,1,3]);end
                
                Imag =  Imag(XPixClip+1:end,YPixClip+1:end,end:-1:1); % clip & flip in z  
                sz =    size(Imag);

                switch modality
                    case 'mus'
                        switch filetype
                            case 'nifti'
                                I =     10.^(Imag/10); 
%                                 
%                                 I = double(I);
% 
%                                 lambda = 1100; %1000:50:1500; %mm: 1480 
%                                 parfor nslice = 1:MZL
% %                                     nslice
%                                     I(:,:,nslice) = denoise_TV_gamma_mm3(sqrt(I(:,:,nslice)),lambda,1,1).^2;
%                                 end
%                                 
%                                 I = single(I);

                                data =  zeros(sz(1),sz(2),MZL,'like',Imag);
                                for z=1:MZL
                                    data(:,:,z) = I(:,:,z)./(sum(I(:,:,z+1:end),3))/2/0.0025;
                                end
                            case 'mat'
                                data = Imag(:,:,1:MZL); 
                        end
                        %mus_3d = zeros(sz, 'like', Imag);
                        %mus_3d(:,:,1:MZL) = data;
                        %for z=MZL+1:sz(3)
                        %    mus_3d(:,:,z) = I(:,:,z)./(sum(I(:,:,z+1:end),3))/2/0.0025;
                        %end
                        %mus_2d = mean(mus_3d(:,:,5:end),3);
                        %fout = [outdir filesep '../mus/' filesep modality '_slice' sprintf('%03i',sliceid_out(ss)) '.mat'];
                        %save(fout, 'mus_2d')
                    case 'dBI'
                        data = Imag(:,:,1:MZL); 
                end
                
                     
                Mosaic(row,columns,:) = Mosaic(row,columns,:)+data.*RampOrig3;
                Masque(row,columns) = Masque(row,columns)+RampOrig;        
            end
            M=M+Mosaic;
            Ma=Ma+Masque;   
        end
    end

    Ma = repmat(Ma, [1,1,MZL]);
    
    MosaicFinal=M./(Ma);
    MosaicFinal(isnan(MosaicFinal))=0;
    MosaicFinal(isinf(MosaicFinal))=0;
    
    if 1==1 % approve this flag when 'MosaicFinal' is smaller than 10 G 
        disp(' - Saving mosaic...');
        fout = [outdir filesep modality '_slice' sprintf('%03i',sliceid_out) '.mat'];
        fprintf(' %s \n',fout);
        aatic=tic;
        save(fout, 'MosaicFinal', 'modality', '-v7.3');
        %parsave(fout, MosaicFinal);
        aatoc=toc(aatic);
        fprintf(' - %.1f s to save\n',aatoc);
    end 
    
%     if strcmpi(modality, 'mus')
%         fout = [outdir filesep '../StitchingFiji' filesep modality '_slice' sprintf('%03i',sliceid_out(ss)) '.mat'];
%         mus_2d = mean(MosaicFinal(:,:,5:end),3);
%         mus_2d = rot90(mus_2d,-1);
%         save(fout, 'mus_2d');
%     end


%reset(gpuDevice())
%%                    
%%%% FOR DOWNSAMPLING
disp(' - Downsampling to 20 micron iso (xy)...');

sx = sxyz(1); %10um -> 20um
sy = sxyz(2); %10um -> 20um
sz = 1; %2.5um stay the same
ds_xyz = [sx, sy, sz];

[Id] = do_downsample(MosaicFinal, ds_xyz);

%         view3D(Id);
fout = [outdir filesep modality '_slice' sprintf('%03i',sliceid_out) '_mean-xy_20um-iso.mat'];
%save(fout, 'Id', 'modality', '-v7.3');
parsave(fout, Id);


toc0=toc(tic0);
disp(['Elapsed time to stitch ' num2str(length(sliceid_out)) ' slices: ']);
disp(['      ' num2str(toc0) ' seconds']);

end


%%
s = 1;
poolobj = parpool(10);
for s = 2:5
    filename = sprintf('/autofs/cluster/octdata2/users/Hui/20220201_SP002/process_run1//mus/mm_despeckle_dBI180/dBI_slice%03i.mat',s);
    load(filename,'MosaicFinal');
%% mm despeckle
I =     10.^(MosaicFinal/10); 

I = double(I);

lambda = 1100; %1000:50:1500; %mm: 1480 

MZL_1 = 50;
parfor nslice = 1:MZL_1% MZL
    nslice
    I(:,:,nslice) = denoise_TV_gamma_mm3(sqrt(I(:,:,nslice)),lambda,1,1).^2;
end

I = single(I);

%% dBI format & saving
tic
dBI_mm = 10*log10(I(:,:,1:MZL_1));
fout = [outdir filesep 'dBI_slice' sprintf('%03i',s) '_mmdespeckle.mat'];
fprintf('saving dBI_mm\n')
save(fout, 'dBI_mm', '-v7.3');
toc
%% mus calculation

tic
sz = size(I);
fprintf('calculating mus\n')
mus =  zeros(sz(1),sz(2),MZL_1,'like',I);
for z=1:MZL_1
    mus(:,:,z) = I(:,:,z)./(sum(I(:,:,z+1:end),3))/2/0.0025;
end

toc

%% save mus
tic
fout = [outdir filesep 'mus_slice' sprintf('%03i',s) '_mmdespeckle.mat'];
fprintf('saving mus\n')
save(fout, 'mus', '-v7.3');
toc
end


%%
list_z = [1:5:46, 50];
list_lambda = [1500:50:2000];
name_out = 'dBI_mm_despeckle_montage.x_z.y_lambda.normalized_w2';
t.decay = squeeze(mean(I3_og_dBI,[1 2]));
t.weight = t.decay(1)./t.decay;

clear time_matrix I_matrix_mm
for iz = 1:length(list_z)
    z = list_z(iz)
    I = I3_og_dBI(:,:,z);
    I = I*t.weight(z);

    I_delog = 10.^(I./10);
    for il = 1:length(list_lambda)
        lambda = list_lambda(il)
        tic
        I_mm = denoise_TV_gamma_mm3(sqrt(I_delog),lambda,1,1).^2; fprintf('.\n')
        time_matrix(il,iz)=toc;
    
        I_matrix_mm{il,iz} = I_mm;

    end
end
% figure; imagesc(10*log10([I_in_mm; I_out_mm_1500; I_out_mm_1700; I_out_mm_1800 ;I_out_mm_1900 ;I_out_mm_2000]),[50 90]); colormap(gray); axis equal
% caxis([60 80])

I = 10*log10(cell2mat(I_matrix_mm));


name = ['/autofs/cluster/octdata2/users/Hui/20220201_SP002/process_run1/dBI/mm-despeckle/20220712_lambda_z_crop/' name_out];
save([name '.mat'],'I','list_z','list_lambda','time_matrix')
imwrite(mat2gray(I,[60 80]),[name '.tiff'],'compression','none');

%% time for one z with lambda (parfor) ...

poolobj = parpool(20);
lambda = 1750
I =     10.^(MosaicFinal/10); 

I = double(I);
tic
parfor nslice = 1:15% MZL
    nslice
    I(:,:,nslice) = denoise_TV_gamma_mm3(sqrt(I(:,:,nslice)),lambda,1,1).^2;
end
toc
I = single(I);

end