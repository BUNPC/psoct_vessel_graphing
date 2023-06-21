function [mosaic_xy] = stitch_xy(mparams, scan, params, modality)
%% Stitch the x-y tiles for each slice.
% This script will determine the x-y coordinates for a single slice of a
% single run. Then, it will apply these coordinates to each subsequent
% slice in the run.
% INPUTS:
%   mparams (struct): mosaic parameters (input dir., output dir., input
%                   file type, slice index, etc.)
%   scan (struct): parameters from scan sequence
%   params (struct): (slice ID, tile ID, gray range, agar threshold)
%   modality (string): mus or dBI
% OUTPUTS:
%   mosaic_xy (matrix): X-Y stitching

%% Initalization params
% Initialize output struct
mosaic_xy = struct;

% slice indices for run
sliceidx    = mparams.sliceidx;

% mosaic parameters
fprintf(' - Loading Experiment file...\n %s\n', mparams.Exp);
S = whos('-file',mparams.Exp);
if any(contains({S(:).name},'Experiment_Fiji'))
    idx = find(contains({S(:).name},'Experiment_Fiji'));
elseif any(contains({S(:).name},'Experiment'))
    idx = find(contains({S(:).name},'Experiment'));
end
load(mparams.Exp,S(idx).name);
Experiment  = eval(S(idx).name);
fprintf(' - %s is loaded ...\n %s\n', S(idx).name);

% input/output directories
indir       = mparams.indir;
outdir      = mparams.outdir;
if ~exist(outdir,'dir')
    mkdir(outdir)
end
filetype    = mparams.InFileType; % 'nifti';
fprintf(' - Input directory = \n%s\n',indir{:});
fprintf(' - Output directory = \n%s\n',outdir);

% Set the x,y clipping (depends on system)
switch scan.System
    case 'Octopus'
        XPixClip    = params.YPixClip; %18;
        YPixClip    = params.XPixClip; %0;
    case 'Telesto'
        XPixClip    = params.XPixClip; %18;
        YPixClip    = params.YPixClip; %0;
end
fprintf('XPixClip = %i\nYPixClip = %i\nSystem is %s\n\n',...
    XPixClip, YPixClip,scan.System);

% Set other parameters
MZL     = 40; %mparams.MZL;
fprintf('save %ipx of orignal %ipx\n', MZL, scan.CropDepth);

%%% 
NbPix = Experiment.NbPix;
X       = Experiment.X_Mean;           Y       = Experiment.Y_Mean;
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

%% Stitching
tabrow=size(Experiment.MapIndex_Tot,1);
tabcol=size(Experiment.MapIndex_Tot,2);

% Loop slices
for s = 1:size(sliceidx,2)
    tic0=tic;
            
    sliceid_in  = sliceidx(1,s); 
    sliceid_out = sliceidx(2,s);
    sliceid_run = sliceidx(3,s); indir_curr = indir{sliceid_run};
    
    MapIndex = Experiment.MapIndex_Tot_offset+Experiment.First_Tile -1;
    fprintf('\nStarting %s mosaic slice %d from run %d slice %d ...\n',...
        modality,sliceid_out, sliceid_run, sliceid_in);

    M  = zeros(MXL,MYL,MZL,'single'); %flipped xy for telesto
    Ma = zeros(MXL,MYL,'single');
    Mosaic = zeros(MXL,MYL,MZL,'single');
    Masque = zeros(MXL,MYL,'single');

    % Loop through columns
    for ii=1:tabcol
        fprintf('col %d / %d\n',ii,tabcol)

        % Loop through rows
        for jj=1:tabrow
            if MapIndex(jj,ii)>0 && ~isnan(X(jj,ii))
                columns =Y(jj,ii):Y(jj,ii)+sizecol-1;
                row     =X(jj,ii):X(jj,ii)+sizerow-1;
                currtile = (sliceid_in-1)*Experiment.TilesPerSlice+MapIndex(jj,ii);
                switch filetype
                    case 'mat'
                        ftile = [indir_curr filesep modality '_' sprintf('%03i',currtilel) '.mat'];
                        S = whos('-file',ftile);
                        load(ftile);
                        Imag = eval(S.name);
                    case 'nifti'
                        ftile = [indir_curr filesep 'test_processed_' sprintf('%03i',currtile) '_cropped.nii'];
                        Imag =  readnifti(ftile);
                end  
                
                Imag = single(Imag);
                if strcmpi(scan.System, 'Octopus');Imag = permute(Imag,[2,1,3]);end              
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
    m_xy= M./(Ma);
    m_xy(isnan(m_xy)) = 0;
    m_xy(isinf(m_xy)) = 0;
    
    %%% Add x-y stitched mus to struct (index by slice #)
    mosaic_xy(s).m_xy = m_xy;
    mosaic_xy(s).modality = modality;

    %%% Save m_xy (raw resolution)
    %{
    disp(' - Saving mosaic...');
    fout = [outdir, filesep, modality, '_slice', sprintf('%03i',sliceid_out), '.mat'];
    fprintf(' %s \n',fout);
    aatic=tic;
    save(fout, 'm_xy', 'modality', '-v7.3');
    aatoc=toc(aatic);
    fprintf(' - %.1f s to save\n',aatoc);
    %}

toc0=toc(tic0);
disp(['Elapsed time to stitch ', num2str(length(sliceid_out)), ' slices: ']);
disp(['      ', num2str(toc0), ' seconds']);

end
end