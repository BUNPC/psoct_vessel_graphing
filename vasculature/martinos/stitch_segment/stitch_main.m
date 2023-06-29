%% Step 1 of pre-processing
%{
This script will stitch each slice.
For each data collection run:
    - change the ParameterFile for each "process_run#"
    - execute this script
To Do:
    - Test step 2 output with MAT-file object

This script should have the following input parameters:
- resolution
- slice thickness
- multirun
- needregistration
- slice intake number (start:stop)
- floating slice index (index of slices w/ artifacts)
- mosaic dimensions for each run (e.g. 10*11 tiles)

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

%% CAA6 Frontal - Initial Parameters
%{
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
%}

%% CAA26 occipital - Initial Parameters
% This was completed in one run, so there is only one parameter file.

% Base directory of data
data_dir = {'/autofs/space/omega_001/users/caa/CAA26_Occipital/Process_caa26_occipital/'};
% imaging modality (PS-OCT)
modality = 'mus';
% Number of slices in each run
run = struct();
run(1).slices = 1:158;
% voxel resolution
res = [0.01 0.01 0.0025];
% voxel downsampling
ds_xyz = [1, 1, 4];
% Define whether needs registration or is uniform
need2register = 'isuniform';
% Use filler for variable regrun since interrun registration is unecessary
regrun = 'None';


%% Step 1: x-y stitching

% CAA 26 - occipital
%{
%%% Initialization Parameters
% Scan Parameters
run1_param_dir = '/autofs/space/omega_001/users/caa/CAA26_Occipital/Process_caa26_occipital/Parameters.mat';
load(run1_param_dir, 'Mosaic3D', 'Scan', 'Parameters');
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
% Create MAT-file object of run1_mosaic_xy
M_xy = matfile(fout,'Writable',true);
%}

% CAA6_frontal
%{
%%% Run 1
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
%}

%% Step 2: z stitching
% TODO: determine if I should create a MAT-file Object of m_xyz to pass
% into the Frangi filter script. Unsure if the Mat file would be too large.
%{
[m_xyz] = stitch_z(M_xy, run1_mparams, res, ds_xyz,...
    run, data_dir, case_notes, regrun, modality);

% Save raw output
fname = strcat(modality, '_mean_10um-iso.slice40px.nii');
fout = fullfile(run1_outdir, fname);
res_ds = res .* ds_xyz;
save_mri(m_xyz, fout, res_ds)

% Save smoothed output
fname = strcat(modality, '_mean_10um-iso.slice40px.sigma1.nii');
fout = fullfile(run1_outdir, fname);
save_mri(mgaussian(m_xyz,1), fout, res_ds, 'float')
%}

%% Step 3: Frangi segmentation


%% Step 4: remove pia

%% Step 5: revise Frangi

%%% Step 6: profits
sprintf('$$$$$$')

%%% WIP code for serializing function call
% ParameterFile = {'', '', ''};
% nruns = length(ParameterFile);
% for run = 1:nruns
%     stitch_xy(ParameterFile{run}, modality);
% end

