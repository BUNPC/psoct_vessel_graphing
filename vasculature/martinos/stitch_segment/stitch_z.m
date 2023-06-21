function [m_xyz] =...
    stitch_z(M_xy, mparams, res, ds_xyz, runs, data_dir, case_notes, regrun, modality)
%STITCH_Z stitch the z axis
% Takes the X-Y mosaic as an input and stitches together adjacent layers in
% the z direction. The output is the reconstructed volume (m_xyz).
%
% INPUTS:
%   M_xy (MAT-file obj.): object connected to the MAT-file for x-y mosaic.
%       This enables memory-efficient access. M_xy will contain an XY
%       mosaic matrix for each run. For example:
%           M_xy.run1_m_xy
%           M_xy.run2_m_xy
%           ...
%           M_xy.run[N]_m_xy
%            
%   mparams (struct): mosaic parameters (input dir., output dir., input
%                   file type, slice index, etc.)
%
%   res (array): voxel resolution [x,y,z]
%
%   ds_xyz (array): downsampling scalar for each axis [x,y,z]
%
%   runs (struct): slice indices for each run. Ex: run(#).slices = 1:50
%
%   data_dir (cell): cell of strings for each directory path to multirun
%       ex: {'/autofs/cluster/octdata2/users/Chao/caa/caa_6/frontal/process_run1/dBI'}
%
%   case_notes (string): state whether interrun registration is necessary
%       ex: 'need2register' or 'isuniform'
%
%   regrun (string): directory/filename for interrun registration
%       ex: [directory]/imregcorr_p.mat
%       This is only needed if interrun registration is needed.
%
%   modality (string): imaging modality
%       'dBI' or 'mus'
%
% OUTPUTS:
%   m_xyz (matrix): x-y-z stitched mosaic

%% TODO:
% add logic to iterate over runs in M_xy

%% Generate a slice indexing matrix (3xN)
% We must generate an array of slice indices for each run. These slice
% indices are stored in the variable "Mosaic3D.sliceidx", which is in the
% file: /autofs/cluster/octdata2/users/Chao/caa/caa_6/frontal/process_run1/Parameters.mat'
%
% mparams.sliceidx (3xN matrix)
%       First row = slice number within run
%       second row = overall slice number across all runs
%       third row = array of ones (always)
%
% sliceidx matrix (3xN matrix)
%    - First row = slice number within run
%    - second row = overall slice number across all runs
%    - third row = run index from data_dir. The run number
%                   may not always start at 1, depending on whether the
%                   user chooses to include run 1 or not.

if length(runs) == 1
    %%% 1 run
    sliceidx = mparams.sliceidx(1,:);
else
    %%% Multiple runs
    % Input slice indices
    s_in = [];
    for ii = 1:length(run)
        s_in = [s_in, run(ii).slices];
    end
    
    % Output slice indices
    s_out = 1:length(s_in);
    
    % Run index
    s_run = [];
    for ii = 1:length(run)
        s_run = [s_run, ones(length(run(ii).slices)).*ii];
    end
        
    % Final slice index matrix
    sliceidx = [s_in; s_out; s_run];
    
    % Load registration parameters (if interrun translation needed)
    if strcmp(case_notes, 'need2register')
        load(regrun,'crop_idx','tran_idx')
    end
end

%% Generate z stitching parameters (multiple runs)
zoff = mparams.z_parameters(1);
z_sm = mparams.z_parameters(2);
z_s = mparams.z_parameters(3);
z_sms       = z_sm+z_s; 
z_m         = z_sm-z_s;

% Print skirt information to console
fprintf('zoff    skirt      mainbody    skirt\n')
fprintf('%0.2i      %0.2i         %0.2i          %0.2i\n',zoff, z_s,z_m,z_s)
fprintf('%0.2i + %0.2i = %0.2ipx (%ium) -> slice thickness\n',z_s,z_m,z_sm,z_sm*2.5)

%% Generate WM mask, calc. tissue decay, calc. & blend weights, 
% This section will use the second slice from first run, since sometimes
% there can be issues with the first slice. This section will access the
% MAT-file object M_xy to more efficiently use memory.

% Get fieldnames of M_xy MAT-file Object
fields = fieldnames(M_xy);
% Take first field name, which corresponds to first run, even in the case
% of a single run.
run1_m_xy = M_xy(fields{1});
% Take second slice from first run
m_xy = run1_m_xy(1,2);

% If using scattering coefficient, then calculate tissue decay
if strcmpi(modality,'mus')
    %%% Calculate average intensity projection (AIP) of XY mosaic slice
    mus_aip = mean(smooth3(m_xy), 3);
    thresh = multithresh(mus_aip,2); % output will be [agar-gm gm-wm]
    figure;subplot(1,2,1);imshow(mus_aip,[thresh(2),prctile(mus_aip(:),98)]); title('wm');
    subplot(1,2,2);imshow(mus_aip,[thresh(1),thresh(2)]); title('gm');

    %%% Mask for white matter
    wm_mask = mus_aip>thresh(2);
    mask3d = logical(m_xy(1,1,:)*0 + wm_mask);
    %%% Mask for grey matter
    % gm_mask = mus_aip<thresh(2) & mus_aip>thresh(1);
    % mask3d = logical(Id(1,1,:)*0 + gm_mask);
    
    %%% Calculate tissue decay
    tissue = mean(reshape(m_xy(mask3d),[],size(m_xy,3)), 1);
    clear mask3d
    tissue = tissue(:).';
    % Plot tissue decay
    figure; plot(tissue); title('tissue decay')

    %%% generate weights with uniform intensity level across depth
    tissue = tissue((1:z_sms)+zoff);
    s  = tissue(z_s+1)./tissue(1:z_s);        % Top Overlapping Skirt
    ms = tissue(z_s+1)./tissue(z_s+1:z_sms);  % Non-overlapping Body  &  Bottom Overlapping Skirt
    
    %%% blend weight at the overlap
    switch modality
        case 'dBI'; degree = 1;   
        case 'mus'; degree = 1; % 1  
    end
    w1 = s             .*linspace(0,1,z_s).^degree; w1 = w1(:).';    % Top Overlapping Skirt
    w2 = ms(z_m+1:end) .*linspace(1,0,z_s).^degree; w2 = w2(:).';    % Bottom Overlapping Skirt
    w3 = ms(1:z_m);                                 w3 = w3(:).';

elseif strcmpi(modality,'conv')
    w1 = linspace(0,1,z_s); w1 = w1(:).';    % Top Overlapping Skirt
    w2 = linspace(1,0,z_s); w2 = w2(:).';    % Bottom Overlapping Skirt
    w3 = (1:z_m)*0+1; w3 = w3(:).';
end

%%% Calculate size of XYZ mosaic
nb_slices   = size(sliceidx,2);
[row_pxblock, col_pxblock] = size(m_xy,[1 2]);
tot_z_pix   = z_sm * nb_slices;
Ma_sz = [row_pxblock, col_pxblock, tot_z_pix];
Ma_sz = Ma_sz ./ ds_xyz;

%%% initialize stitched xyz mosaic
m_xyz = zeros(Ma_sz,'single');
clear m_xy s ms

%%% Set input/output directories
% Print input directory
for ii = 1:length(data_dir)
    fprintf(' - Input directory = \n%s\n',data_dir{ii})
end

% Set output directory
outdir = mparams.outdir;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
fprintf(' - Output directory = \n%s\n',outdir);



%% Begin stitching  
% Idea:
%   - add outer for-loop to iterate over runs. 
%   - create counter for slice index
%   - load slice from run
% 
tic
% Iterate over slices
for s = 1:length(sliceidx)
    s_in = sliceidx(1,s);
    s_out = sliceidx(2,s);
    s_run = sliceidx(3,s);
    fprintf('\tslice %d in run %d \tslice %d in all runs\n',s_in,s_run,s_out);
    indir = data_dir{s_run};

    % TODO: load slice from MAT-file object
    

    fdata = [indir filesep modality '_slice' sprintf('%03i',s_in) '.mat'];
    load(fdata,'MosaicFinal');
    MosaicFinal = MosaicFinal(:, :, zoff + 1:z_sms);
    ddims = size(MosaicFinal);

    switch case_notes
        case 'need2register'
            c1s=crop_idx(i_run,1);
            c1e=crop_idx(i_run,2);
            c2s=crop_idx(i_run,3);
            c2e=crop_idx(i_run,4);
            t1s=tran_idx(i_run,1);
            t1e=(c1e-c1s) + tran_idx(i_run,1); % crop length + translation index
            t2s=tran_idx(i_run,2);
            t2e=(c2e-c2s) + tran_idx(i_run,2);
        case 'isuniform'
            c1s=1;
            c1e=ddims(1);
            c2s=1;
            c2e=ddims(2);
            t1s=1;
            t1e=ddims(1); % crop length + translation index
            t2s=1;
            t2e=ddims(2);
    end

    %%% stitch (cases for first, last, and other slice #s)
    if s_out==1 % first_slice
        W = reshape(repmat([w1*0+1 w3 w2],[ddims(1)*ddims(2),1]),...
            ddims(1),ddims(2),[]);
        zr1 = 1:z_sms / ds_xyz(3);
        zr2 = 1:z_sms + zoff;
        % Downsample
        I_in = do_downsample(MosaicFinal(:,:,zr2).*W, ds_xyz);
        m_xyz(t1s:t1e,t2s:t2e,zr1) = I_in(c1s:c1e,c2s:c2e,:);
    elseif s_out==length(sliceidx) % last_slice
        % Add to the `t` bottom-most pixels in m_xyz
        W = reshape(repmat([w1 w3 ],[ddims(1)*ddims(2),1]),...
            ddims(1),ddims(2),[]);            
        zr1 = (1:z_sm/ds_xyz(3)) + Ma_sz(3)-z_sm/ds_xyz(3);
        zr2 = zoff+(1:z_sm);
        % Downsample
        I_in = do_downsample(MosaicFinal(:,:,zr2).*W, ds_xyz);
        m_xyz(t1s:t1e,t2s:t2e,zr1) = m_xyz(t1s:t1e,t2s:t2e,zr1) + I_in(c1s:c1e,c2s:c2e,:);
    else
        % For slices (first_slice+1) to (last_slice-1)
        W = reshape(repmat([w1 w3 w2],[ddims(1)*ddims(2),1]),...
            ddims(1),ddims(2),[]);
        zr1 = (1:z_sms/ds_xyz(3)) + (s-1)*z_sm/ds_xyz(3);
        zr2 = zoff+(1:z_sms);
        % Downsample
        I_in = do_downsample(MosaicFinal(:,:,zr2).*W, ds_xyz);
        m_xyz(t1s:t1e,t2s:t2e,zr1) = m_xyz(t1s:t1e,t2s:t2e,zr1) + I_in(c1s:c1e,c2s:c2e,:);
    end
    disp('Finished slice');
end

toc
end