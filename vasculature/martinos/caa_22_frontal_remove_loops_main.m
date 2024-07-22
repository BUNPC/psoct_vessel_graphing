%% Remove loops main script
% Purpose: remove the loops from the graph of CAA 17

clear; close all; clc;

%% Add top-level directory of code repository to path
% Start in current directory
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

%% Import the segmentation
% Path to masked segmentation
dpath = '/projectnb/npbssmic/ns/CAA/caa22/frontal/segmentations/';
seg_name = 'caa22-frontal_vessels-masked.mat';
% Import segmentation
seg = load(fullfile(dpath, seg_name));
f = fields(seg);
seg = seg.(f{1});

%% Convert segmentation to graph and save
% Voxel dimensions (microns)
vox_dim = [20, 20, 20];

% Boolean for visualizing the graph debugging plots
viz = false;

% Boolean for removing loops from graph
rmloop_bool = 1;

% Output file path and filename of output (without suffix)
fullpath = '/projectnb/npbssmic/ns/CAA/caa22/frontal/segmentations/';
fname = 'caa22-frontal_vessels-masked';

% Create a graph of the segmentation
seg_graph_init(seg, vox_dim, fullpath, fname, viz, rmloop_bool); 