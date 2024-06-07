%% Main file for calling segmentation functions
% Author: Mack Hyman
% Date Created: June 05, 2024
%
% Detailed Description
%{
This script performs the following:
- Import the CAA vascular and EPVS segmentation
- Convert segmentation to a graph (remove loops)
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

%% Generate skeletons/graphs of cropped volumes
% Use cropped regions of subjects CAA22-occipital and CAA22-frontal.
% Compare metrics (length density, branch density, tortuosity) from the
% skeleton of the EPVS and vascular segmentation
%{
%%% CAA 22 Occipital (cropped 350x350xZ voxels)
caa22_path = fullfile(dpath,'caa22/occipital/segmentations/');
seg = fullfile(caa22_path,'caa22-occipital_vessels-masked_crop_350.tif');
seg = logical(TIFF2MAT(seg));
epvs = fullfile(caa22_path,'epvs_crop_350.tif');
epvs = logical(TIFF2MAT(epvs));

% VESSELS & EPVS: create graph and skeleton (with loops)
viz = false;
rm_loop = false;
vessel_fname = 'caa22_occipital_vessels_crop_350';
seg_graph_init(seg,[20,20,20],caa22_path,vessel_fname,viz,rm_loop);
epvs_fname = 'caa22_occipital_epvs_crop_350';
seg_graph_init(epvs,[20,20,20],caa22_path,epvs_fname,viz,rm_loop);

% VESSELS & EPVS: create graph (remove loops)
% note that the saved file with have "rm_loops" in the name, which will
% distinguish it from the former file
rm_loop = true;
seg_graph_init(seg,[20,20,20],caa22_path,vessel_fname,viz,rm_loop);
seg_graph_init(epvs,[20,20,20],caa22_path,epvs_fname,viz,rm_loop);
%}

%%% CAA 22 Frontal (cropped 350x350xZ voxels)
caa22_path = fullfile(dpath,'caa22/frontal/segmentations/');
seg = fullfile(caa22_path,'caa22-frontal_vessels-masked_crop_350.tif');
seg = logical(TIFF2MAT(seg));
epvs = fullfile(caa22_path,'epvs_min15vox_crop_350.tif');
epvs = logical(TIFF2MAT(epvs));

% VESSELS & EPVS: create graph and skeleton (with loops)
viz = false;
rm_loop = false;
vessel_fname = 'caa22-frontal_vessels-masked_crop_350';
epvs_fname = 'caa22_frontal_epvs_min15vox_crop_350';
% seg_graph_init(seg,[20,20,20],caa22_path,vessel_fname,viz,rm_loop);
% seg_graph_init(epvs,[20,20,20],caa22_path,epvs_fname,viz,rm_loop);

% VESSELS & EPVS: create graph (remove loops)
% note that the saved file with have "rm_loops" in the name, which will
% distinguish it from the former file
rm_loop = true;
sprintf('Initializing Graph for CAA 22 frontal vessels (no loops)\n')
seg_graph_init(seg,[20,20,20],caa22_path,vessel_fname,viz,rm_loop);
sprintf('Finished Graph for CAA 22 frontal vessels (no loops)\n')
sprintf('Initializing Graph for CAA 22 frontal EPVS (no loops)\n')
seg_graph_init(epvs,[20,20,20],caa22_path,epvs_fname,viz,rm_loop);
sprintf('Finished Graph for CAA 22 frontal EPVS (no loops)\n')