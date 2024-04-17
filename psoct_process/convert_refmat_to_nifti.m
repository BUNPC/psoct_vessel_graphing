%% Convert mask.mat and ref.mat to NIFTI format
% Each subdirectory ".../[subjectID]/dist_corrected/volume/ref contains the
% c-scans and mask c-scans in the format ref.mat and mask.mat. The mask
% must be corrected in freeview to remove false positives/negatives. This
% script will iterate through these files for each subject and create a
% NIFTI for each.

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

%% Initialize folders and file names
% Subfolder containing data
subdir = '/dist_corrected/volume/ref';
% Filename to parse (this will be the same for each subject)
fbase = 'ref';
%%% Complete subject ID list for Ann_Mckee_samples_10T
subid = {'AD_20832', 'AD_20969',...
         'AD_21354', 'AD_21424',...
         'CTE_6489','CTE_6912',...
         'CTE_7019','CTE_8572','CTE_7126',...
         'NC_6839','NC_6974','NC_8653',...
         'NC_21499','NC_301181'};

% Filenames
refname = 'ref.mat';
maskname = 'mask.mat';
% Specifics for NC_6047
%{
subid = {'NC_6047'};
refname = 'refc.tif';
%}

% for ii = 1:length(subid)
parfor (ii = 1:length(subid), NSLOTS)
    %% Load "ref.mat" and "mask.mat" from the subfolder
    % Define entire filepath 
    refpath = fullfile(dpath, subid{ii}, subdir, refname);
    maskpath = fullfile(dpath, subid{ii}, subdir, maskname);
    % Load ref and mask
    if contains(refname, 'ref.mat')
        ref = load(refpath);
        mask = load(maskpath);
        ref = ref.vol;
        mask = mask.vol;
    else
        % Import ref
        ref = TIFF2MAT(refpath);
        mask = load(maskpath);
        mask = mask.tiss_mat;
    end

    %% Initialize filename and header information
    % Names of output files
    refout = fullfile(dpath, subid{ii}, subdir, 'ref.nii');
    maskout = fullfile(dpath, subid{ii}, subdir, 'mask.nii');

    % Resolution of volumes
    res = [0.012, 0.012, 0.015];

    % Datatypes
    ref_dtype = 'double';
    mask_dtype = 'uchar';

    % Permute flag (1 = swap x/y, 0 = don't)
    permuteflag = 1;
    
    %% Save each as a nifti
    save_mri(ref, refout, res, ref_dtype, permuteflag)
    save_mri(mask, maskout, res, mask_dtype, permuteflag)
end