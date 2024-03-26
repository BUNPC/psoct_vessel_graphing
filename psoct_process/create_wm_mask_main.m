%% Main script to create the white matter mask
% This main script calls the function create_wm_mask(mus,th), which
% thresholds the scattering map (mus) by the threshold (th). This function
% uses a GUI for manually modfiying the mask to remove any false positives
% and fill in false negatives.

%% Create a mask for the white matter
% Use the scattering map to segment just the white matter
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
vol_name = 'ref';

%% Set maximum number of threads
% Retrieve the number of available cores
n_cores = str2num(getenv('NSLOTS'));
% Set the maximum number of threads equal to the number of cores
maxNumCompThreads(n_cores);

% Instantiate parallel pool (if needed)
% poolobj = gcp('nocreate');
% if isempty(poolobj)
%     % No parallel pool exists. Initialize a parallel pool.
%     pc = parcluster('local');
%     % Setup directory for logging
%     pc.JobStorageLocation = getenv('TMPDIR');
%     % Start running the parallel pool
%     parpool(pc, n_cores);
% end

%% Initialize subject ID lists
%%% All subjects to analyze
subid = {'AD_10382', 'AD_20832', 'AD_20969',...
         'AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912',...
         'CTE_7019', 'CTE_7126',...
         'NC_6839',  'NC_6974', 'NC_8653',...
         'NC_21499', 'NC_301181'};

%%% Number of scattering maps (mus_#.tif) per volume
% Note, the scattering maps for CTE_8572 must be generated
nmus = [17,15,22,15,20,...
        22,10,18,18,...
        16,20,15,15,24];
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
    %%% Debugging information
    fprintf('\n---------Starting Subject %s---------\n',subid{ii})

    %%% Import non-normalized volume
    fprintf('importing non-normalized volume\n')
    fpath = fullfile(dpath, subid{ii}, subdir, strcat(vol_name,'.mat'));
    vol = load(fpath); 
    vol = vol.vol;

    %%% Retrieve last slice # and N slices per physical stack
    sid = char(subid{ii});
    % Retrieve the last slice #
    zmax = regstruct.(sid).zmax;
    
    %%% Initialize matrices for storing masks
    wm_mat = zeros(size(vol));
    gm_mat = zeros(size(vol));
    tiss_mat = zeros(size(vol));

    %% Create white mask for each physical slice   
    % Counter for physical slice
    mus_cnt = 1;
    for j=1:nstack:zmax
        %%% Scattering map frame
        fprintf('Importing mus frame %i \n',mus_cnt)
        % Create string of scattering map
        mus_str = strcat('mus',num2str(mus_cnt),'.tif');
        fpath = fullfile(dpath, sid, scatdir, mus_str);
        % Rescale the scattering map to uint8 for registration
        mus = uint8(rescale(TIFF2MAT(fpath),0,255));
        % Iterate the mus slice count
        mus_cnt = mus_cnt + 1;
        
        %%% OCT intensity frame
        % The scattering map has extra padding of zeros along the bottom
        % and right borders, compared to the OCT frame. This section of
        % code will remove these zeros from the borders.
        fprintf('Importing OCT frame %i \n',j)
        oct = vol(:,:,j);
        % Find difference in matrix size with scattering map
        x = size(mus,1) - size(oct,1);
        y = size(mus,2) - size(oct,2);
        % Truncate the zero paddings from x and y
        mus = mus(1:end-x, 1:end-y);

        %%% Segment the white matter from the mus
        % Overlay the oct and mus
        figure; h = imshow(oct);
        set(h, 'AlphaData',mus);
        % Call function to apply white matter mask
        th = 110;
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

    % TODO: 
    % - bitwise multiply masks and refined mask
    % - save the mask to the subject's folder
    % - save the WM/GM to the subject's folder


end


















