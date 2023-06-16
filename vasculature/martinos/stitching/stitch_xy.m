function [MosaicFinal] = stitch_xy(ParameterFile, modality)
%% Stitch the x-y tiles for each slice.
% This script will determine the x-y coordinates for a single slice of a
% single run. Then, it will apply these coordinates to each subsequent
% slice in the run.
% INPUTS:
%   ParameterFile (string): 
%   modality (string): 
% OUTPUTS:
%   MosaicFinal (matrix): X-Y stitching

%%% Load Mosaic3D parameters variable
load(ParameterFile, 'Mosaic3D', 'Scan', 'Parameters');

%%% SETTING SLICE INPUT 
sliceidx    = Mosaic3D.sliceidx;

%%% SETTING MOSAIC PARAMETERS 
fprintf(' - Loading Experiment file...\n %s\n', Mosaic3D.Exp);
S = whos('-file',Mosaic3D.Exp);
if any(contains({S(:).name},'Experiment_Fiji'))
    idx = find(contains({S(:).name},'Experiment_Fiji'));
elseif any(contains({S(:).name},'Experiment'))
    idx = find(contains({S(:).name},'Experiment'));
end
load(Mosaic3D.Exp,S(idx).name);
Experiment  = eval(S(idx).name);
fprintf(' - %s is loaded ...\n %s\n', S(idx).name);

%%% SETTING DIRECTORIES 
indir       = Mosaic3D.indir;
outdir      = Mosaic3D.outdir;
if ~exist(outdir,'dir')
    mkdir(outdir)
end
filetype    = Mosaic3D.InFileType; % 'nifti';

fprintf(' - Input directory = \n%s\n',indir{:});
fprintf(' - Output directory = \n%s\n',outdir);

% Set the x,y clipping (depends on system)
switch Scan.System
    case 'Octopus'
        XPixClip    = Parameters.YPixClip; %18;
        YPixClip    = Parameters.XPixClip; %0;
    case 'Telesto'
        XPixClip    = Parameters.XPixClip; %18;
        YPixClip    = Parameters.YPixClip; %0;
end
fprintf('XPixClip = %i\nYPixClip = %i\nSystem is %s\n\n', XPixClip, YPixClip,Scan.System);

% Set other parameters
MZL     = 40; %Mosaic3D.MZL;
fprintf('save %ipx of orignal %ipx\n', MZL, Scan.CropDepth);

%%% 
NbPix = Experiment.NbPix;

X       = Experiment.X_Mean;           Y       = Experiment.Y_Mean;         % Y = fliplr(Y);
X       = X-min(X(:))+1;               Y       = Y-min(Y(:))+1;
sizerow = NbPix-XPixClip;              sizecol = NbPix-YPixClip;
MXL     = max(X(:))+sizerow-1;         MYL     = max(Y(:))+sizecol-1;

%%% GENERATE RAMP
ramp_x = sizerow-round(median(diff(Experiment.X_Mean,[],1),'all','omitnan'));
ramp_y = sizecol-round(median(diff(Experiment.Y_Mean,[],2),'all','omitnan'));

ramp_xv       = 0:ramp_x-1;            ramp_yv       = 0:ramp_y-1;
x             = ones(1,sizerow);       y             = ones(1,sizecol);
x(1+ramp_xv)  = mat2gray(ramp_xv);     y(1+ramp_yv)  = mat2gray(ramp_yv);
x(end-ramp_xv)= mat2gray(ramp_xv);     y(end-ramp_yv)= mat2gray(ramp_yv);

RampOrig=x.'*y;
RampOrig3 = repmat(RampOrig, [1,1,MZL]);

%% 
tabrow=size(Experiment.MapIndex_Tot,1);
tabcol=size(Experiment.MapIndex_Tot,2);

% Loop slices
for s = 1:size(sliceidx,2)
    tic0=tic;
            
    sliceid_in  = sliceidx(1,s); 
    sliceid_out = sliceidx(2,s);
    sliceid_run = sliceidx(3,s); indir_curr = indir{sliceid_run};
    
    MapIndex = Experiment.MapIndex_Tot_offset+Experiment.First_Tile -1;
    fprintf('\nStarting %s mosaic slice %d from run %d slice %d ...\n',modality,sliceid_out, sliceid_run, sliceid_in);

    M  = zeros(MXL,MYL,MZL,'single'); %flipped xy for telesto
    Ma = zeros(MXL,MYL,'single');
    Mosaic = zeros(MXL,MYL,MZL,'single');
    Masque = zeros(MXL,MYL,'single');
    
    for ii=1:tabcol
        fprintf('col %d / %d\n',ii,tabcol)
        for jj=1:tabrow
            if MapIndex(jj,ii)>0 && ~isnan(X(jj,ii))
                columns =Y(jj,ii):Y(jj,ii)+sizecol-1;
                row     =X(jj,ii):X(jj,ii)+sizerow-1;
                currtile = (sliceid_in-1)*Experiment.TilesPerSlice+MapIndex(jj,ii);
                switch filetype
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
                if strcmpi(Scan.System, 'Octopus');Imag = permute(Imag,[2,1,3]);end              
                Imag =  Imag(XPixClip+1:end,YPixClip+1:end,end:-1:1); % clip & flip in z  
                sz =    size(Imag);

                switch modality
                    case 'mus'
                        switch filetype
                            case 'nifti'
                                I =     10.^(Imag/10); 

                                data =  zeros(sz(1),sz(2),MZL,'like',Imag);
                                for z=1:MZL
                                    data(:,:,z) = I(:,:,z)./(sum(I(:,:,z+1:end),3))/2/0.0025;
                                end
                            case 'mat'
                                data = Imag(:,:,1:MZL); 
                        end
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

    %%% Final stitching
    Ma = repmat(Ma, [1,1,MZL]);
    MosaicFinal=M./(Ma);
    MosaicFinal(isnan(MosaicFinal))=0;
    MosaicFinal(isinf(MosaicFinal))=0;
    
    %%% Save MosaicFinal (raw resolution)
    disp(' - Saving mosaic...');
    fout = [outdir filesep modality '_slice' sprintf('%03i',sliceid_out) '.mat'];
    fprintf(' %s \n',fout);
    aatic=tic;
    save(fout, 'MosaicFinal', 'modality', '-v7.3');
    aatoc=toc(aatic);
    fprintf(' - %.1f s to save\n',aatoc);

    % Save the 2D projection
    if strcmpi(modality, 'mus')
        fout = [outdir filesep '../StitchingFiji' filesep modality '_slice' sprintf('%03i',sliceid_out(ss)) '.mat'];
        mus_2d = mean(MosaicFinal(:,:,5:end),3);
        mus_2d = rot90(mus_2d,-1);
        save(fout, 'mus_2d');
    end

toc0=toc(tic0);
disp(['Elapsed time to stitch ' num2str(length(sliceid_out)) ' slices: ']);
disp(['      ' num2str(toc0) ' seconds']);

end
end