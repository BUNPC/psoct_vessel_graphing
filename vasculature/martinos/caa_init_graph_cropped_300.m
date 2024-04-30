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
topdir = mydir(1:idcs(end-1));
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
fnames = {'caa6-frontal_vessels-masked_crop_300',...
          'caa6-occipital_vessels-masked_crop_300',...
          'caa17_occipital_THRESH-0.5_crop_300',...
          'caa22-frontal_vessels-masked_crop_300',...
          'caa25-frontal_vessels-masked_crop_300',...
          'caa25-occipital_vessels-masked_crop_300',...
          'caa26-frontal_vessels-masked_crop_300',...
          'caa26-occipital_vessels-masked_crop_300'};
% filename extension
ext = '.tif';

%%% Complete subject ID list for Ann_Mckee_samples_10T
subids = {'caa6/frontal/', 'caa6/occipital/',...
         'caa17/occipital/',...
         'caa22/frontal/',...
         'caa25/frontal/', 'caa25/occipital/',...
         'caa26/frontal/', 'caa26/occipital/'};

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
seg = TIFF2MAT(filename);

%% Convert to graph and save
% Voxel dimensions (microns)
vox_dim = [20, 20, 20];

% Boolean for visualizing the graph debugging plots
viz = false;

% Boolean for removing loops from graph
rmloop_bool = 1;

%%% Create a graph of the segmentation
seg_graph_init(seg, vox_dim, fullpath, fname, viz, rmloop_bool); 