%% Main file for calling segmentation functions
% Author: Mack Hyman
% Date Created: April 24, 2024
%
% Detailed Description
%{
This script performs the following:
- Import the CAA segmentation (nifti file)
- Convert segmentation to a logical
- Save the segmentation
- Convert segmentation to a graph
%}
clear; clc; close all;

%% Add top-level directory of code repository to path
% This allows Matlab to find the functions in the project folders
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

%% Set maximum number of cores
% Retrieve the number of available cores
n_cores = str2num(getenv('NSLOTS'));
% Set the maximum number of threads equal to the number of cores
maxNumCompThreads(n_cores);

%% Initialize data path for linux or personal machine (debugging)
% Path to top-level directory
dpath = '/projectnb/npbssmic/ns/CAA/';
% Subfolder containing data
subdir = 'segmentations/';
% Filename to parse (this will be the same for each subject)
fnames = {'caa6-frontal_vessels-masked',...
          'caa6-occipital_vessels-masked',...
          'caa17-occipital_vessels-masked',...
          'caa22-frontal_vessels-masked',...
          'caa25-frontal_vessels-masked',...
          'caa25-occipital_vessels-masked',...
          'caa26-frontal_vessels-masked',...
          'caa26-occipital_vessels-masked'};
% filename extension
ext = '.nii';

%%% Complete subject ID list for Ann_Mckee_samples_10T
subids = {'caa_6/frontal/', 'caa_6/occipital/',...
         'caa_17/occipital/',...
         'caa_22/frontal/',...
         'caa_25/frontal/', 'caa_25/occipital/',...
         'caa_26/frontal/', 'caa_26/occipital/'};

%%% Reassign subid and sigma based on job array counter
% Retrieve SGE_TASK_ID from system (job array index)
batch_idx = getenv('SGE_TASK_ID');
% Convert from ASCII to double
batch_idx = str2double(batch_idx);
% Retrieve corresponding row from sub_sigma
subid = subids{batch_idx};
fname = fnames{batch_idx};

%% Load segmentation (NII) and convert to MAT
% Define entire filepath 
fullpath = fullfile(dpath, subid, subdir);
filename = strcat(fullpath, strcat(fname, ext));
% Convert .tif to .MAT
seg = MRIread(filename, 0, 0);
seg = logical(seg.vol);

%% Save Segmentation
% Save masked segmentation as a sparse matrix in .MAT
fout = strcat(fullfile(fullpath, fname), '.mat');
save(fout, 'seg', '-v7.3');

%% Convert to graph and save

% Voxel dimensions
vox_dim = [20, 20, 20];

% Boolean for visualizing the graph debugging plots
viz = false;

% Boolean for removing loops from graph
rmloop_bool = 0;

%%% Create a graph of the segmentation
seg_graph_init(seg, vox_dim, fullpath, fname, viz, rmloop_bool); 