function Mosaic2D_Telesto(ParameterFile, modality, nWorker, method)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIAL SETTING 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(ParameterFile); %Parameters Mosaic2D Scan

%%% if using FOVAdjust (hard coded to false)
FovAdjustFlag = 0;

%%% if using resize_img 
resize_img  = Mosaic2D.imresize;

%%% if method specified 
method_select = 0;
if nargin > 3; method_select = 1; end 
%%% display numbers of workers 
% fprintf('--nWorker is %i--\n', nWorker);

%%% display modality 
fprintf('--Modality is %s--\n', modality);

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
sliceidx    = Mosaic2D.sliceidx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING MOSAIC PARAMETERS 

fprintf(' - Loading Experiment file...\n %s\n', Mosaic2D.Exp);
S = whos('-file',Mosaic2D.Exp);
if any(contains({S(:).name},'Experiment_Fiji'))
    idx = find(contains({S(:).name},'Experiment_Fiji'));    isFiji = true;
elseif any(contains({S(:).name},'Experiment'))
    idx = find(contains({S(:).name},'Experiment'));         isFiji = false;
end
load(Mosaic2D.Exp,S(idx).name);
Experiment  = eval(S(idx).name);
fprintf(' - %s is loaded ...\n %s\n', S(idx).name);

if method_select==1
    Experiment.X_Mean = Experiment.(method).X_Mean;
    Experiment.Y_Mean = Experiment.(method).Y_Mean;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING DIRECTORIES 
indir       = Mosaic2D.indir;
outdir      = Mosaic2D.outdir; 
if method_select==1; outdir = [outdir filesep method]; end
if ~exist(outdir,'dir'); mkdir(outdir); end
filetype    = Mosaic2D.InFileType; % 'nifti';

%fprintf(' - Input directory = \n%s\n',indir{:});
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING IMAGE GRAYSCALE RANGE
GrayRange = [];
savetiff  = 1;
switch lower(modality)
    case 'aip';         if isfield(Parameters,'AipGrayRange'); GrayRange = Parameters.AipGrayRange; else; noscale = 0;end
    case 'mip';         if isfield(Parameters,'MipGrayRange'); GrayRange = Parameters.MipGrayRange; else; noscale = 0;end
    case 'retardance';  if isfield(Parameters,'RetGrayRange'); GrayRange = Parameters.RetGrayRange; else; noscale = 0;end
    case 'mus';         if isfield(Parameters,'musGrayRange'); GrayRange = Parameters.musGrayRange; else; noscale = 0;end
    otherwise;          if isfield(Parameters,[modality 'GrayRange']); GrayRange = Parameters.([modality 'GrayRange']); else; noscale = 0;end %disp(' - unrecognized modality!');
end
if exist('noscale', 'var')
    savetiff = 0; 
    warning(' %s Grayscale is not found. Will only output MAT file. ', modality);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 
NbPix = Experiment.NbPix;

X       = Experiment.X_Mean;           Y       = Experiment.Y_Mean;         % Y = fliplr(Y);
X       = X-min(X(:))+1;               Y       = Y-min(Y(:))+1;
sizerow = NbPix-XPixClip;              sizecol = NbPix-YPixClip;
MXL     = max(X(:))+sizerow-1;         MYL     = max(Y(:))+sizecol-1;

if strcmpi(modality,'orientation');MZL=4;else;MZL=1;end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set mosaic params
%MapIndex = Experiment.MapIndex_Tot_offset+ Experiment.First_Tile - 1;


tabrow=size(Experiment.MapIndex_Tot,1);
tabcol=size(Experiment.MapIndex_Tot,2);

sample = double(1); % single(1);
sample_whos = whos('sample');
fprintf(' - Processing in %s\n',sample_whos.class);
fprintf(' - Saving in %s\n','single');

%poolobj = parpool(nWorker);
%parfor s = 1:size(sliceidx,2)
for s = 1:size(sliceidx,2)

            
    sliceid_in  = sliceidx(1,s); 
    sliceid_out = sliceidx(2,s);
    sliceid_run = sliceidx(3,s); 
    indir_curr = indir{sliceid_run}; % this code can accommdate datas from multiple runs
%    indir_curr = indir;

%     if isFiji
%         MapIndex = squeeze(Experiment.MapIndex_Tot(:,:,sliceid_out));
%     else
        MapIndex = Experiment.MapIndex_Tot_offset+Experiment.First_Tile -1;
%     end
    % Start mosaic
    tic
    fprintf('\nStarting %s mosaic slice %d from run %d slice %d ...\n',modality,sliceid_out, sliceid_run, sliceid_in);

   

    
    M   =zeros(MXL,MYL,MZL,'like',sample);
    Ma  =zeros(MXL,MYL,'like',sample);

    for ii=1:tabcol
    %     fprintf('Column %d\n',ii);
        for jj=1:tabrow
            Mosaic  =zeros(MXL,MYL,MZL,'like',sample);
            Masque  =zeros(MXL,MYL,'like',sample);

    %         fprintf('row %d\n',jj);
            if MapIndex(jj,ii)>0 && ~isnan(X(jj,ii))
                columns  = Y(jj,ii):Y(jj,ii)+sizecol-1;
                row      = X(jj,ii):X(jj,ii)+sizerow-1;

%                 if isFiji;  currtile = MapIndex(jj,ii);
%                 else;       currtile = (sliceid_in-1)*Experiment.TilesPerSlice+MapIndex(jj,ii); 
%                 end
                currtile = (sliceid_in-1)*Experiment.TilesPerSlice+MapIndex(jj,ii);
                
                % load/generate I (2D image)
                switch filetype
                    case 'nifti'
                        if strcmpi(modality,'mus')
                            Imag=[indir_curr filesep 'test_processed_' sprintf('%03i',currtile) '_cropped.nii'];
                            I3 = readnifti(Imag); 
                            if ~isa(I3, 'single'); I3 = single(I3);end
                            if canUseGPU(); I3 =  gpuArray(I3);end
                            if strcmpi(Scan.System, 'Octopus');I3 = permute(I3,[2,1,3]);end

                            I3 = I3(XPixClip+1:end,YPixClip+1:end,:);
                            sz = size(I3);
                            I3 = I3(:,:,end:-1:1);
                            I3 = 10.^(I3/10); 

                            data =  zeros(sz,'like',I3);
                            for z=1:sz(3)-1
                                data(:,:,z) = I3(:,:,z)./(sum(I3(:,:,z+1:end),3))/2/0.0025;
                            end
                            I = squeeze(mean(data(:,:,5:70),3)); %figure;plot(squeeze(mean(data,[1 2])))
                        else % orientation, aip, mip, retardance
                            Imag=[indir_curr filesep 'test_processed_' sprintf('%03i',currtile) '_' lower(modality) '.nii'];
                            I = niftiread(Imag);
                            %if ~isa(I, 'single'); I = single(I);end
                            if strcmpi(Scan.System, 'Octopus');I = I.';end
                            I = I(XPixClip+1:end,YPixClip+1:end);
                        end
                    case 'mat'
                        Imag=[indir_curr filesep modality '_' sprintf('%03i',currtile) '.mat'];
                            S = whos('-file',Imag);
                            load(Imag);
                            I = eval(S.name);
                            if strcmpi(Scan.System, 'Octopus');I = I.';end
                            I = I(XPixClip+1:end,YPixClip+1:end);
                end

                % blend I into Mosaic
                if strcmpi(modality,'orientation')
                    %I = orientationSign * I + orientationOffset;
                    SLO_OC=zeros([size(I) 4],'double');
                    SLO_OC(:,:,1)=cos(I/180*pi).*cos(I/180*pi);
                    SLO_OC(:,:,2)=cos(I/180*pi).*sin(I/180*pi);
                    SLO_OC(:,:,3)=SLO_OC(:,:,2);
                    SLO_OC(:,:,4)=sin(I/180*pi).*sin(I/180*pi);

                    Mosaic(row,columns,:)=SLO_OC.*repmat(RampOrig,1,1,4);
                    Masque(row,columns) = RampOrig;
                else % mus, aip, mip, retardance
                    Mosaic(row,columns) = I.*RampOrig;
                    Masque(row,columns) = RampOrig;
                end
            end
            M=M+Mosaic;
            Ma=Ma+Masque;
        end
    end

    switch lower(modality)
        case 'mip';         modalstr = 'MIP';
        case 'aip';         modalstr = 'AIP';
        case 'retardance';  modalstr = 'Retardance';
        case 'orientation'; modalstr = 'Orientation';
        case 'mus';         modalstr = 'mus';
        otherwise;          modalstr = modality;%error('foo:bar','Unknown image modality!');
    end


    if strcmpi(modality, 'orientation')
        fprintf('Starting orientation angles eigen decomp...\n');
        tic
        for st =  1:size(M,1)
    %             fprintf('...\n');
            for tt = 1:size(M,2)
                % calculate tensor, and do eigen decomposition
                a = [M(st,tt,1) M(st,tt,3);M(st,tt,2) M(st,tt,4)];
                [V, D] = eig(a);
                Data2_1(st,tt) =  atan(V(2,2)/V(1,2));
            end
        end
        toc
        
        %%% orientation angles eigen decomposition
        M=M./(Ma);
        M = double(M);
        a_x = reshape(M,[size(M,1)*size(M,2),2,2]);
        [V,~]= eig2(a_x);
        Data_ = atan(V(:,2,2)./V(:,1,2));
        Data_ = reshape(Data_, size(MosaicFinal,[1 2]));
        Data  = Data_/pi*180;
        O     = Data;

        index1=O<-90;
        index2=O>90;
        O(index1)=O(index1)+180;
        O(index2)=O(index2)-180;
        O=rot90(O,-1);
        MosaicFinal=O;
        

        
        %%% save MosaicFinal
        if ~isa(MosaicFinal, 'single'); MosaicFinal = single(MosaicFinal);end
        fprintf('Saving Orientation mosaic .mat...\n'); 
        parsave([outdir filesep  modalstr '_slice' sprintf('%03i',sliceid_out) '.mat'], MosaicFinal); % ,'-v7.3');

        %%% masking orientation
        RetSliceTiff = [outdir filesep 'Retardance_slice' sprintf('%03i',sliceid_out) '.tiff'];
        AipSliceTiff = [outdir filesep 'AIP_slice' sprintf('%03i',sliceid_out) '.tiff'];
        data1 = imread(RetSliceTiff);
        data2 = imread(AipSliceTiff);
        data4 = wiener2(data2,[5 5]);

        map3D=ones([size(O,[1,2]) 3]);
        map3D2=ones([size(O,[1,2]) 3]);
        O_normalized = (O+90)/180;

        %%% Orientation1
        I=mat2gray(double(data1))/(1-0.4); 
    %    I=(double(data1)-80)/60; % -lowerthreshold) / DynamicRange
        I(I<0)=0;
        I(I>1)=1;
        I(data4<=20)=0;
        map3D(:,:,1)=O_normalized;
        map3D(:,:,3)=I;
        maprgb=hsv2rgb(map3D);
    %         figure; imshow(maprgb); title('Orientation 1');
        fprintf('Saving Orientation1 mosaic .tiff...\n');
        imwrite(maprgb,[outdir filesep  modalstr '1_slice' sprintf('%03i',sliceid_out) '.tiff'],'compression','none');
        if resize_img
            fprintf('Resizing image...\n');
            inimg = imread([outdir filesep  modalstr '1_slice' sprintf('%03i',sliceid_out) '.tiff']);
            imrsz = imresize(inimg,0.3,'nearest');
            % figure; imshow(imrsz); title('resized image');
            imwrite(imrsz,[outdir filesep  'resized,nearest-' modalstr '1_slice' sprintf('%03i',sliceid_out) '.tiff'],'compression','none');
        end

        %%% Orientation2
        I=(-mat2gray(double(data2))+1-0.2)/((1-0.2)*0.5); % -() +1 invert; -0.2 clip graymatter; (1-0.2) then dynamic range; 0.5 bring shadow to midtone
    %    I=(double(data2)+20)/220;
        I(I<0)=0;
        I(I>1)=1;
        I(data4<=20)=0;
        map3D2(:,:,1)=O_normalized;
        map3D2(:,:,3)=I;
        maprgb2=hsv2rgb(map3D2);
    %     figure; imshow(maprgb2); title('Orientation 2');
        fprintf('Saving Orientation2 mosaic .tiff...\n');
        imwrite(maprgb2,[outdir filesep modalstr '2_slice' sprintf('%03i',sliceid_out) '.tiff'],'compression','none');
        if resize_img
            fprintf('Resizing Orientation2 image...\n');
            inimg = imread([outdir filesep modalstr '2_slice' sprintf('%03i',sliceid_out) '.tiff']);
            imrsz = imresize(inimg,0.3,'nearest');
            % figure; imshow(imrsz); title('resized image');
            imwrite(imrsz,[outdir filesep  'resized,nearest-' modalstr '2_slice' sprintf('%03i',sliceid_out) '.tiff'],'compression','none');
        end

    else
        % If not Orientation, save tiff and mat for slice
        MosaicFinal=rot90(M./(Ma),-1);
        MosaicFinal(isnan(MosaicFinal))=0;


        if savetiff == 1
            fprintf('Saving .tiff mosaic...\n');
            foutimg = [outdir filesep modalstr '_slice' sprintf('%03i',sliceid_out) '.tiff'];
            imwrite(mat2gray(MosaicFinal,GrayRange),foutimg,'compression','none');
            
            if resize_img == 1
                fprintf('Resizing image...\n');
                inimg = imread(foutimg);
                imrsz = imresize(inimg,0.5,'nearest');
                % figure; imshow(imrsz); title('resized image');
                foutimg2 = [outdir filesep 'resized-' modalstr '_slice' sprintf('%03i',sliceid_out) '.jpg'];
                imwrite(imrsz,foutimg2);
            end
        end
        
        if ~isa(MosaicFinal, 'single'); MosaicFinal = single(MosaicFinal);end
        fprintf('Saving .mat mosaic...\n');
        foutmat = [outdir filesep modalstr '_slice' sprintf('%03i',sliceid_out) '.mat'];
        parsave(foutmat, MosaicFinal);

    end


    toc
    fprintf('\n   ---   Finished with %s slice %d   ---   \n',modalstr,sliceid_out);
end

delete(gcp('nocreate'));
end







% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

