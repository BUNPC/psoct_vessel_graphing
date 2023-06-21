%% Step 1 of pre-processing
%{
This script will stitch each slice.
For each data collection run:
    - change the ParameterFile for each "process_run#"
    - execute this script
To Do:
    - parameterize function call for each run with for-loop

% Manually verify that there are no overlapping slies.
% 1) Navigate to the directory "StitchingFiji"
% 2) Using Fiji, open the following images:
%     - slice(end) from run001
%     - slice001 in run002
% 3) Verify these images are from different slices
% 4) Perform this validation for all subsequent runs.

% Processed Directory
% Note that this will vary for each subject, depending on who ran the first
% processing. Here is one example of the subfolder:
% /autofs/cluster/octdata2/users/Chao/caa/[subID]/processed/[YYYYMMDD]/[YYYYMMDD]_[run#]_processed/
% /autofs/cluster/octdata2/users/Chao/caa/[subID]/proces_run#
% look at the parameters.mat filepath structures

%}
clear; clc;

%% Add top-level directory to path
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Truncate path to reach top-level directory (psoct_vessel_graphing)
topdir = mydir(1:idcs(end-2));
addpath(genpath(topdir));

%% Initial Parameters
% Base directory of data
base_dir = '/autofs/cluster/octdata2/users/Chao/caa/caa_6/frontal';

% Subfolder containing processed runs
proc_dir = fullfile(base_dir, 'process_run');

% imaging modality (PS-OCT)
modality = 'mus';

% Number of slices in each run
run = struct();
run(1).slices = 1:50;
run(2).slices = 1:14;
run(3).slices = 1:14;

% voxel resolution
res = [0.01 0.01 0.0025];
% voxel downsampling
ds_xyz = [1, 1, 4];

% Multirun process directory
multirun_p.ProcDir = {'/autofs/cluster/octdata2/users/Chao/caa/caa_6/frontal/process_run1',...
    '/autofs/cluster/octdata2/users/Chao/caa/caa_6/frontal/process_run2',...
    '/autofs/cluster/octdata2/users/Chao/caa/caa_6/frontal/process_run3'};

%% Step 1: x-y stitching
% TODO: use a for-loop to iterate over runs

%%% Run 1
%{
% Initialization Parameters
run1_params = '/autofs/cluster/octdata2/users/Chao/caa/caa_6/frontal/process_run1/Parameters.mat';
load(run1_params, 'Mosaic3D', 'Scan', 'Parameters');
run1_outdir = Mosaic3D.outdir;
run1_mparams = Mosaic3D;
run1_scan = Scan;
run1_params = Parameters;
% Stitch X-Y
[run1_mosaic_xy] = stitch_xy(run1_mparams, run1_scan, run1_params, modality);
% Save output for debugging;
fname = strcat(modality, '_xy', '.mat');
fout = fullfile(run1_outdir, fname);
save(fout, 'run1_mosaic_xy', 'modality', '-v7.3');
% Convert to MAT-file to access portions of data without loading entirety
run1_mxy = matfile(fout,'Writable',true);
%}

%%% Run 2
% Initialization Parameters
run2_params = '/autofs/cluster/octdata2/users/Chao/caa/caa_6/frontal/process_run2/Parameters.mat';
load(run2_params, 'Mosaic3D', 'Scan', 'Parameters');
run2_outdir = Mosaic3D.outdir;
run2_mparams = Mosaic3D;
run2_scan = Scan;
run2_params = Parameters;
% Stitch X-Y
[run2_mosaic_xy] = stitch_xy(run2_mparams, run2_scan, run2_params, modality);
% Save output for debugging;
fname = strcat(modality, '_xy', '.mat');
fout = fullfile(run2_outdir, fname);
save(fout, 'run2_mosaic_xy', 'modality', '-v7.3');
% Convert to MAT-file to access portions of data without loading entirety
run2_mxy = matfile(fout,'Writable',true);

%%% Run 3
% Initialization Parameters
run3_params = '/autofs/cluster/octdata2/users/Chao/caa/caa_6/frontal/process_run3/Parameters.mat';
load(run3_params, 'Mosaic3D', 'Scan', 'Parameters');
run3_outdir = Mosaic3D.outdir;
run3_mparams = Mosaic3D;
run3_scan = Scan;
run3_params = Parameters;
% Stitch X-Y
[run3_mosaic_xy] = stitch_xy(run3_mparams, run3_scan, run3_params, modality);
% Save output for debugging;
fname = strcat(modality, '_xy', '.mat');
fout = fullfile(run3_outdir, fname);
save(fout, 'run3_mosaic_xy', 'modality', '-v7.3');
% Convert to MAT-file to access portions of data without loading entirety
run3_mxy = matfile(fout,'Writable',true);


%% Step 2: z stitching
% Inputs:
%   - MosaicFinal (from step 1)
% create separate function for "save_mri". Remove this function from end of
% script.

% Number of slices per run (stitch in z)

[run1_mosaic_xyz] = stitch_z(run1_mosaic_xy);
% Pass output to step 3

%%% Step 3: Frangi segmentation

%%% Step 4: remove pia

%%% Step 5: revise Frangi

%%% Step 6: profits
sprintf('$$$$$$')

%%% WIP code for serializing function call
% ParameterFile = {'', '', ''};
% nruns = length(ParameterFile);
% for run = 1:nruns
%     stitch_xy(ParameterFile{run}, modality);
% end

