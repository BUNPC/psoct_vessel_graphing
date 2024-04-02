%% Create mask for tissue and pia boundary
% This script will find the edge of the tissue for each slice, create a
% continuous polygon to outline the tissue, and then erode the boundary of
% the tissue. This will be repeated for each slice. This mask will then be
% applied to the segmentation to remove false positives along the border of
% the tissue.

%% Initialize environment
clear; clc; close all;

% Debug argument (true = display figures for each slice)
debug = true;

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
%     poolobj = gcp('nocreate');
%     if isempty(poolobj)
% 	    myCluster=parcluster('local');
% 	    % Ensure multiple parpool jobs don't overwrite other's temp files
% 	    myCluster.JobStorageLocation = getenv('TMPDIR');
% 	    poolobj = parpool(myCluster, NSLOTS);
%     end

end

% Subfolder containing data
subdir = '/dist_corrected/volume/ref/';
% Filename to parse (this will be the same for each subject)
fname = 'refc.tif';
%%% Complete subject ID list for Ann_Mckee_samples_10T
subid = {'AD_10382', 'AD_20832', 'AD_20969',...
         'AD_21354', 'AD_21424',...
         'CTE_6489','CTE_6912',...
         'CTE_7019','CTE_8572','CTE_7126',...
         'NC_6839','NC_6974',...
         'NC_8653','NC_21499','NC_301181'};
subid = {'NC_6047'};

%% Iterate subjects
for ii = 1:length(subid)
% parfor (ii = 1:length(subid), NSLOTS)
    %%% Debugging information
    fspec = 'Started subject %s\n';
    fprintf(fspec, subid{ii})

    %% Determine number of ref#.mat files
    % Define file path to ref#.mat files
    fullpath = fullfile(dpath, subid{ii}, subdir);
    % Create struct from directory contents
    list = dir(fullpath);
    % Extract names
    names = {list.name};
    % Create regular expression TODO: update this to have 1-2 digits
    exp = 'ref\d+.mat';
    % Find strings matching regexp
    refnames = regexp(names, exp, 'match');
    % Remove empty elements (did not match)
    refnames = refnames(~cellfun('isempty',refnames));
    % Counter for number of ref files
    nref = length(refnames);
    
    %% Initialize matrix to store mask for entire volume
    % Load the volume "ref.mat"
    fname = fullfile(fullpath, 'ref.tif');
    vol = TIFF2MAT(fname);
    % Create uint8 matrix to store mask
    mask = uint8(zeros(size(vol)));
    
    %% Iterate over ref#.mat files and create mask for each stack
    % Initialize counter for location in stack
    z0 = 1;
    % Iterate through ref#.mat stacks and create mask for each
    for j=1:nref
        %%% Load z-stack (TIF) and convert to MAT
        % Update fpath for next ref stack
        fname = strcat('ref',num2str(j),'.mat');
        fname = fullfile(fullpath, fname);
        ref = load(fname);
        if isfield(ref,'ref')
            ref = ref.ref;
        elseif isfield(ref, 'Ref')
            ref = ref.Ref;
        else
            error('Unexpected structure.')
        end
           
        %%% Create mask for each volume
        % Function to create mask for each slice in stack
        [ref_mask, t] = create_mask_v4(ref, debug);

        %%% Add stack mask to volume mask stack
        zf = z0+size(ref,3)-1;
        mask(:,:,z0:zf) = ref_mask;
        z0 = zf + 1;
    end

    %% Overlay tissue volume with active contour mask
    % Call function to apply mask
    volm = apply_mask(vol, mask);
    
    %% Save mask and masked volume    
    % Create output filenames
    mask_out = fullfile(fullpath, strcat('mask'));
    volm_out = fullfile(fullpath, strcat('ref_masked'));
    % Save output as .TIF
    segmat2tif(mask, strcat(mask_out, '.tif'));
    segmat2tif(volm, strcat(volm_out, '.tif'));
    % Save output as .MAT
    save_ref(mask_out, mask);
    save_ref(volm_out, volm);

    %%% Debugging information
    fspec = 'Finished subject %s\n';
    fprintf(fspec, subid{ii})
end