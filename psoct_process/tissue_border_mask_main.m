%% Create mask for tissue and pia boundary
% This script will find the edge of the tissue for each slice, create a
% continuous polygon to outline the tissue, and then erode the boundary of
% the tissue. This will be repeated for each slice. This mask will then be
% applied to the segmentation to remove false positives along the border of
% the tissue.

clear; clc; close all;

%% Add top-level directory of code repository to path
% This allows Matlab to find the functions in the project folders

% Start in current directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Truncate path to reach top-level directory (psoct_vessel_graphing)
topdir = mydir(1:idcs(end));
addpath(genpath(topdir));

%% Initialize data path for linux or personal machine (debugging)

%%% Local machine
if ispc
    dpath = ['C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures'...
            '\test_data\Ann_Mckee_samples_10T\'];
    NSLOTS = 0;
%%% Computing cluster (SCC)
elseif isunix
    % Path to top-level directory
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
    % Set # threads = # cores for job
    NSLOTS = str2num(getenv('NSLOTS'));
    maxNumCompThreads(NSLOTS);

    % Check to see if we already have a parpool, if not create one with
    % our desired parameters
    poolobj = gcp('nocreate');
    if isempty(poolobj)
	    myCluster=parcluster('local');
	    % Ensure multiple parpool jobs don't overwrite other's temp files
	    myCluster.JobStorageLocation = getenv('TMPDIR');
	    poolobj = parpool(myCluster, NSLOTS);
    end

end

% Subfolder containing data
subdir = '/dist_corrected/volume/';
% Filename to parse (this will be the same for each subject)
fname = 'ref.tif';
%%% Complete subject ID list for Ann_Mckee_samples_10T
subid = {'AD_10382', 'AD_20832', 'AD_20969',...
         'AD_21354', 'AD_21424',...
         'CTE_6489','CTE_6912',...
         'CTE_7019','CTE_8572','CTE_7126',...
         'NC_6047', 'NC_6839',...
         'NC_6974', 'NC_7597',...
         'NC_8095', 'NC_8653',...
         'NC_21499','NC_301181'};
subid = {'AD_10382'};

%% Iterate subjects
for ii = 1:length(subid)
    %% Load raw volume (TIF) and convert to MAT
    % Define entire filepath 
    fullpath = fullfile(dpath, subid{ii}, subdir);
    % Define file path
    fpath = fullfile(fullpath, fname);
    % Import the tissue volume TIF stack
    vol = TIFF2MAT(fpath);
       
    %% Create mask for each layer in volume
    % TODO: include inner for-loop to iterate over each "ref#.mat" stack.
    % Using the consolidated "ref.mat" stack is unpredictable since the
    % number of images per physical slice varies.
%     for j=1:nrefs      

    % Debug argument (true = display figures for each slice)
    debug = false;
    % Function to create mask for each slice in stack
    [mask, t] = create_mask_v3(vol, debug, NSLOTS);

    %% Overlay tissue volume with active contour mask
    % Call function to apply mask
    volm = apply_mask(vol, mask);
    
    %% Save mask and masked volume    
    % Create output filenames
    mask_out = fullfile(fullpath, strcat('mask_v3'));
    volm_out = fullfile(fullpath, strcat('ref_4ds_masked_v3'));
    % Save output as .TIF
    segmat2tif(mask, strcat(mask_out, '.tif'));
    segmat2tif(volm, strcat(volm_out, '.tif'));
    % Save output as .MAT
    save(mask_out, 'mask', '-v7.3');
    save(volm_out, 'volm', '-v7.3');
    
end