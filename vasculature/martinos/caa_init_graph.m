%% Convert CAA segmentation to graph and remove loops
% Author: Mack Hyman
% Date Created: April 24, 2024
%
% Detailed Description
%{
This script performs the following:
- Import the CAA segmentation (.nii)
- Import the respective tissue mask (.nii)
- Apply tissue mask
    - save masked vasculature as .MAT
- Convert segmentation to a graph
    - remove loops and save with loops removed
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

%% Initialize data path for segmentations
% Excluding CAA22-frontal because it was already run separately

% Path to top-level directory
dpath = '/projectnb/npbssmic/ns/CAA/';
% Subfolder for subject_ID/region, which contains the mask
subids = {'caa6/frontal/', 'caa6/occipital/',...
         'caa17/occipital/',...
         'caa22/occipital/',...
         'caa25/frontal/', 'caa25/occipital/',...
         'caa26/frontal/', 'caa26/occipital/'};
% Filename of mask
mask_names = {'caa6_frontal_mask_edited.nii',...
          'caa6_occipital_mask_edited.nii',...
          'caa17_occipital_mask_5x-dilate.nii',...
          'caa22_occipital_mask_edited.nii',...
          'caa25_frontal_mask_edited.nii',...
          'caa25_occipital_mask_edited.nii',...
          'caa26_frontal_mask_edited.nii',...
          'caa26_occipital_mask.nii'};
% Subfolder containing segmentation
vasc_dir = 'segmentations/';
% Filename of vasculature segmentation
vasc_names = {'caa6_frontal_THRESH-0.5.nii',...
    'caa6_occipital_THRESH-0.5.nii',...
    'caa17_occipital_THRESH-0.5.nii',...
    'caa22-occipital_prediction-r2-BINARY.nii',...
    'caa25_frontal_THRESH-0.5.nii',...
    'caa25_occipital_THRESH-0.5.nii',...
    'caa26_frontal_THRESH-0.5.nii',...
    'caa26_occipital_THRESH-0.5.nii'};
% Output filename
output_names = {'caa6_frontal_vessels_edited_masked',...
          'caa6_occipital_vessels_edited_masked',...
          'caa17_occipital_vessels_edited_masked',...
          'caa22_occipital_vessels_edited_masked',...
          'caa25_frontal_vessels_edited_masked',...
          'caa25_occipital_vessels_edited_masked',...
          'caa26_frontal_vessels_edited_masked',...
          'caa26_occipital_vessels_edited_masked'};

%%% Reassign subid and sigma based on job array counter
% Retrieve SGE_TASK_ID from system (job array index)
batch_idx = getenv('SGE_TASK_ID');
if strcmp(batch_idx,'undefined')
    batch_idx = 3;
else
    % Convert from ASCII to double
    batch_idx = str2double(batch_idx);
end
% Retrieve corresponding row from sub_sigma
subid = subids{batch_idx};
mask_name = mask_names{batch_idx};
seg_name = vasc_names{batch_idx};

%% Load segmentation and tissue mask
% Segmentation
seg_path = fullfile(dpath, subid, vasc_dir);
filename = fullfile(seg_path, seg_name);
seg = MRIread(filename,0,0);
seg = logical(seg.vol);
% Mask
mask_path = fullfile(dpath, subid);
filename = fullfile(mask_path, mask_name);
mask = MRIread(filename,0,0);
mask = logical(mask.vol);
% Truncate the dimensions for CAA17/occipital
if batch_idx == 3
    mask = mask(:,1:1549,1:596);
    seg = seg(:,1:1549,1:596);
end
% Apply mask to segmentation
seg_masked = logical(mask .* seg);
% Save the masked segmentation prior to loop removal
fname_out = output_names{batch_idx};
fout = fullfile(seg_path, append(fname_out,'.mat'));
save(fout,'seg_masked','-v7.3');

%% Convert to graph and save
% Voxel dimensions (microns)
vox_dim = [20, 20, 20];

% Boolean for visualizing the graph debugging plots
viz = false;

% Boolean for removing loops from graph
rmloop_bool = 1;

%%% Create a graph of the segmentation
seg_graph_init(seg_masked, vox_dim, seg_path, fname_out, viz, rmloop_bool); 