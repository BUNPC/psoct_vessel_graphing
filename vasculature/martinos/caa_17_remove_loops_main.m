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

%% Import the graph and segmentation
% Path to graph
dpath = '/projectnb/npbssmic/ns/CAA/caa_17/occipital/caa17-occipital_masked-filtered/seg_masked_graph.mat';
% Load graph
graph = load(dpath);
graph = graph.graph;
nodes = graph.nodes;
edges = graph.edges;
vox = graph.vox;

% Path to segmentation
seg_path = '/projectnb/npbssmic/ns/CAA/caa_17/occipital/caa17-occipital_masked-filtered/seg_masked.mat';
vol = load(seg_path);
% Import segmentation
vol = vol.seg_masked;

%% Input parameters for Remove Loops Parallel function 
% Move to mean minimum voxel intensity
v_min = 0.99;
% Down sample search radius
delta = 4;
% # iterations for mv2mean function in for-loop iteration in rm_loops
mv_iter = 1;
% Boolean for visualizing debugging graphs
viz = false;

% Call the function to remove loops
[nodes_rm, edges_rm] =...
    rm_loops_parallel(nodes, edges, vol, delta, v_min, mv_iter, viz);

% Save the output to the same directory
graph.nodes = nodes_rm;
graph.edges = edges_rm;
fout = '/projectnb/npbssmic/ns/CAA/caa_17/occipital/caa17-occipital_masked-filtered/seg_masked_graph_rmloops.mat';
save(fout, 'graph', '-v7.3');