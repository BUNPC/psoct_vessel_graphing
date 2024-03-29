%% Main script to create the white matter mask
% This main script calls the function create_wm_mask(mus,th), which
% thresholds the scattering map (mus) by the threshold (th). This function
% uses a GUI for manually modfiying the mask to remove any false positives
% and fill in false negatives.
%
% Remaining work:
%{
- bitwise multiply masks and refined mask
- save the mask to the subject's folder
- save the WM/GM to the subject's folder
%}
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

% Debugging variable. True to view figures.
viz = true;

%% Create white mask for subjects

for ii = 1:length(subid)
% for ii = 9:length(subid)
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
    wm_mat = false(size(ref,1),size(ref,2),zmax);
    gm_mat = false(size(ref,1),size(ref,2),zmax);
    tiss_mat = false(size(ref,1),size(ref,2),zmax);

    %% Create white mask for each physical slice   
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
        figure; h = imshow(oct);
        set(h, 'AlphaData',mus);
        % Set the mus threshold for white matter
        th = threshold(ii);
        % Call function to apply white matter mask
        [wm_mask, gm_mask, tiss_mask] = create_wm_mask(mus,oct,th);
        
        %%% Save the masks
        % Replicate the masks "nstack" times in z dimension. This will
        % replicate the mask for each depth in the physical slice.
        wm_mask = repmat(wm_mask,1,1,nstack);
        gm_mask = repmat(gm_mask,1,1,nstack);
        tiss_mask = repmat(tiss_mask,1,1,nstack);
        % Add the masks to the matrix
        wm_mat(:,:,j:j+nstack-1) = wm_mask;
        gm_mat(:,:,j:j+nstack-1) = gm_mask;
        tiss_mat(:,:,j:j+nstack-1) = tiss_mask;
    end
    
    %% Combine tissue mask and WM mask
    t_mask = fullfile(dpath,sid,subdir,'maskec.mat');
    if isfile(t_mask)
        % Load tissue mask and truncate to match dimensions of mask
        t_mask = load(t_mask);
        t_mask = t_mask.mask;
        t_mask = logical(t_mask(:,:,1:zmax));
        % Bitwise multiply WM mask and refined tissue mask
        wm_mask = bsxfun(@times,t_mask,cast(wm_mat,'like',t_mask));
        % Bitwise multiply ~WM mask and refined tissue mask
        gm_mask = bsxfun(@times,t_mask,cast(~wm_mat,'like',t_mask));     
    end
   
    %%% Save the output
    mdir = fullfile(dpath,sid,subdir,'/masks/');
    if ~isfolder(mdir)
        % Create the folder to store the masks
        mkdir(mdir);
    end
    % Save outputs as .MAT
    fprintf('Saving outputs as .MAT\n')
    wm_output = fullfile(mdir,'wm_mask.mat');
    gm_output = fullfile(mdir,'gm_mask.mat');
    t_output = fullfile(mdir,'t_mask.mat');
    save(wm_output,"wm_mat",'-v7.3');
    save(gm_output,"gm_mat",'-v7.3');
    save(t_output,"tiss_mat",'-v7.3');
    % Save outputs as .TIF
    fprintf('Saving outputs as .TIF\n')
    wm_output = fullfile(mdir,'wm_mask.tif');
    gm_output = fullfile(mdir,'gm_mask.tif');
    t_output = fullfile(mdir,'t_mask.tif');
    segmat2tif(im2uint8(wm_mask),wm_output)
    segmat2tif(im2uint8(gm_mask),gm_output);
    segmat2tif(im2uint8(t_mask),t_output);
    % Save outputs as .NII
    res = [0.012, 0.012, 0.015];
    ref_dtype = 'uchar';
    permuteflag = 1;
    fprintf('Saving outputs as .NII\n')
    save_mri(wm_mask, fullfile(mdir,'wm_mask.nii'), res, ref_dtype, permuteflag);
    save_mri(gm_mask, fullfile(mdir,'gm_mask.nii'), res, ref_dtype, permuteflag);
    save_mri(t_mask, fullfile(mdir,'t_mask.nii'), res, ref_dtype, permuteflag);

end


















