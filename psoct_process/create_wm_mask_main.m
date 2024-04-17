%% Main script to create the white matter mask
% This main script calls the function create_wm_mask(mus,th), which
% thresholds the scattering map (mus) by the threshold (th). This function
% uses a GUI for manually modfiying the mask to remove any false positives
% and fill in false negatives.
%
% This script outputs three masks: mask_tiss, mask_gm_, and mask_wm, which
% are masks for the tissue, white matter, and gray matter, respectively. 
% The tissue mask is generated in a different script, and it is then
% manually refined in Freeview.
% The gray matter (GM) mask is created by taking the bitwise operation
% between the tissue mask and the inverse of the white matter (WM) mask.

clear; clc; close all;

%% Add top-level directory of code repository to path
% Print current working directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Remove the two sub folders to reach parent
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));

%% Initialize directories and filenames
%%% Directories and file name
% Upper level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T';    
% Subfolder containing ref files
subdir = '/dist_corrected/volume/ref';
% Scattering subdirectory
scatdir = 'fitting_4x/';
% non-normalized PSOCT volume
vol_name = 'ref1';

%% Set maximum number of threads
% Retrieve the number of available cores
n_cores = str2double(getenv('NSLOTS'));
% Set the maximum number of threads equal to the number of cores
maxNumCompThreads(n_cores);

%% Debugging Variables
% Debugging figures visualization
viz = false;
% Boolean to load in the cleaned tissue masks
tmask_bool = true;

%% Initialize subject ID lists
%%% All subjects to analyze
subid = {'AD_10382', 'AD_20832', 'AD_20969',...
         'AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912',...
         'CTE_7019', 'CTE_7126','CTE_8572',...
         'NC_6839',  'NC_6974', 'NC_8653',...
         'NC_21499', 'NC_301181'};

%%% Number of scattering maps (mus_#.tif) per volume
% Note, the scattering maps for CTE_8572 must be generated
nmus = [17,15,22,15,20,...
        22,10,18,18,21,...
        16,20,15,15,24];
% white matter intensity threshold for each subject
threshold = [110,110,75,75,110,...
             110,110,110,110,80,...
             70,110,110,70,110];
% Number of frames in z-dimension per physical slice
nstack = 11;
% Number of frames in z-dimension per volume. This is dependent on the
% number of mus maps that were processed for each tissue volume.
zmax = nmus .* nstack;

% Create dictionary to store last image in stack
regstruct = struct;
for ii = 1:length(subid)
    sub = char(subid(ii));
    regstruct.(sub).zmax = zmax(ii);
end

%% Create white mask for subjects
for ii = 1:length(subid)
    %%% Debugging information
    fprintf('\n---------Starting Subject %s---------\n',subid{ii})

    %%% Import single slice of non-normalized volume
    fprintf('importing non-normalized volume\n')
    fpath = fullfile(dpath, subid{ii}, subdir, strcat(vol_name,'.mat'));
    ref = load(fpath); 
    try
        ref = ref.Ref;
    catch
        ref = ref.ref;
    end

    %%% Retrieve last slice # and N slices per physical stack
    sid = char(subid{ii});
    % Retrieve the last slice #
    zmax = regstruct.(sid).zmax;
    
    %%% Initialize matrices for storing masks
    mask_wm = false(size(ref,1),size(ref,2),zmax);
    if ~tmask_bool
        mask_gm = false(size(ref,1),size(ref,2),zmax);
        mask_tiss = false(size(ref,1),size(ref,2),zmax);
    end
    
    %%% Initialize directory for storing masks
    mdir = fullfile(dpath,sid,subdir,'/masks/');

    %% Create white matter mask for each physical slice   
    % Load one depth of the OCT. This is to measure the size of the OCT
    % matrix, which varies from the size of the mus.
    oct = ref(:,:,1);
    % Counter for physical slice
    mus_cnt = 1;
    for j=1:nstack:zmax
        %%% Scattering map frame
        fprintf('Importing mus frame %i \n',mus_cnt)
        % Create string of scattering map
        if strcmp(sid,'NC_301181')
            mus_str = strcat('mus',num2str(mus_cnt),'_ds4x.tif');
            fpath = fullfile(dpath, sid, scatdir, mus_str);
        elseif strcmp(sid,'CTE_8572')
            mus_str = strcat('bfg',num2str(mus_cnt),'_ds4x.tif');
            fpath = fullfile(dpath, sid, 'fitting_bfg_10x', mus_str);
        else
            mus_str = strcat('mus',num2str(mus_cnt),'.tif');
            fpath = fullfile(dpath, sid, scatdir, mus_str);
        end
        % Rescale the scattering map to uint8 for registration
        mus = uint8(rescale(TIFF2MAT(fpath),0,255));
        % Iterate the mus slice count
        mus_cnt = mus_cnt + 1;
        
        %%% OCT intensity frame
        % The scattering map has extra padding of zeros along the bottom
        % and right borders, compared to the OCT frame. This section of
        % code will remove these zeros from the borders.
        x = size(mus,1) - size(oct,1);
        y = size(mus,2) - size(oct,2);
        % Truncate the zero paddings from x and y
        mus = mus(1:end-x, 1:end-y);

        %%% Segment the white matter from the mus
        % Overlay the oct and mus
        if viz
            figure; h = imshow(oct);
            set(h, 'AlphaData',mus);
        end
        % Set the mus threshold for white matter
        th = threshold(ii);
        % Call function to apply white matter mask
        if tmask_bool
            [mask_wm_tmp, ~, ~] = create_wm_mask(mus,oct,th,viz);
        else
            [mask_wm_tmp, mask_gm_tmp, mask_tiss_tmp] =...
                create_wm_mask(mus,oct,th,viz);
        end
        
        %%% Add masks to matrix
        % Replicate the masks "nstack" times in z dimension. This will
        % replicate the mask for each depth in the physical slice.
        mask_wm_tmp = repmat(mask_wm_tmp,1,1,nstack);
        % Add the masks to the matrix
        mask_wm(:,:,j:j+nstack-1) = mask_wm_tmp;
        
        % Repeat for gray matter and tissue masks (commented out because
        % they are generated later on by bitwise operation with the tissue
        % mask.
        if ~tmask_bool
            mask_gm_tmp = repmat(mask_gm_tmp,1,1,nstack);
            mask_tiss_tmp = repmat(mask_tiss_tmp,1,1,nstack);
            mask_gm(:,:,j:j+nstack-1) = mask_gm_tmp;
            mask_tiss(:,:,j:j+nstack-1) = mask_tiss_tmp;
        end
    end
    
    %% Combine tissue mask and WM mask
    
    %%% Check for cleaned white matter mask
    mask_wm_path = fullfile(mdir,'wm_maskec.nii');
    if isfile(mask_wm_path)
        mask_wm = MRIread(mask_wm_path,0,1);
        mask_wm = mask_wm.vol;
        mask_wm = logical(mask_wm(:,:,1:zmax));
    end

    %%% Replace automated tissue mask with the cleaned tissue mask
    if tmask_bool
        mask_tiss = fullfile(dpath,sid,subdir,'maskec.nii');
        % Load tissue mask and truncate to match dimensions of mask
        mask_tiss = MRIread(mask_tiss,0,1);
        mask_tiss = mask_tiss.vol;
        mask_tiss = logical(mask_tiss(:,:,1:zmax));
        % Bitwise operation to generate GM
        mask_gm = bsxfun(@times,mask_tiss ,cast(~mask_wm,'like',mask_tiss));
        mask_gm = logical(mask_gm);
    end
   
    %%% Save the output
    % Create the folder to store the masks
    if ~isfolder(mdir)
        mkdir(mdir);
    end
    
    % Save outputs as .MAT
    fprintf('Saving outputs as .MAT\n')
    wm_output = fullfile(mdir,'mask_wm.mat');
    gm_output = fullfile(mdir,'mask_gm.mat');
    t_output = fullfile(mdir,'mask_tiss.mat');
    save(wm_output,"mask_wm",'-v7.3');
    save(gm_output,"mask_gm",'-v7.3');
    save(t_output,"mask_tiss",'-v7.3');
    
    % Save outputs as .TIF
    fprintf('Saving outputs as .TIF\n')
    wm_output = fullfile(mdir,'mask_wm.tif');
    gm_output = fullfile(mdir,'mask_gm.tif');
    t_output = fullfile(mdir,'mask_tiss.tif');
    segmat2tif(im2uint8(mask_wm),wm_output)
    segmat2tif(im2uint8(mask_gm),gm_output);
    segmat2tif(im2uint8(mask_tiss),t_output);
    
    % Save outputs as .NII
    res = [0.012, 0.012, 0.015];
    ref_dtype = 'uchar';
    permuteflag = 1;
    fprintf('Saving outputs as .NII\n')
    save_mri(mask_wm, fullfile(mdir,'mask_wm.nii'), res, ref_dtype, permuteflag);
    save_mri(mask_gm, fullfile(mdir,'mask_gm.nii'), res, ref_dtype, permuteflag);
    save_mri(mask_tiss, fullfile(mdir,'mask_tiss.nii'), res, ref_dtype, permuteflag);

end