%% Test rm_loops_parallel function
%{
Purpose: The function "rm_loops" operates on each loop in series. This
results in significant processing time when datasets contain large amounts
of loops. This test script is for developing and testing a function to
remove loops in parallel.
%}
clear; close all; clc;

%% Flag for visualization while debugging
visual = false;

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

%% Initialize data paths for dataset with loops
%%% AD_10382
% Top-level directories
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
subid = 'NC_21499';
subdir = '/dist_corrected/volume/';
sigdir = '/combined_segs/gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/';
vdata = 'ref_4ds_norm_inv_refined_masked.tif';
% Data with many nested loops (probability threshold = 0.21)
seg_name = 'seg_refined_masked.tif';
gdata = 'seg_refined_masked_graph_data.mat';

%%% Output Data filenames
skel_out = strcat(gdata(1:end-4),'_loops_rm.tif');
skel_out = char(fullfile(dpath, subid, subdir, sigdir, skel_out));
gdata_out = strcat(gdata(1:end-4),'_loops_rm.mat');
%% Load PSOCT graph, volume, segmentation

%%% Load Graph
Data = load(fullfile(dpath, subid, subdir, sigdir, gdata), 'Data');
Data = Data.Data;
nodes = Data.Graph.nodes;
edges = Data.Graph.edges;

%%% Load volumetric information and set threshold
% Import volume
vol = TIFF2MAT(fullfile(dpath, subid, subdir, vdata));

%% Overlay graph and segmentation
if visual
    %%% Load segmentation stack
    seg = TIFF2MAT(fullfile(dpath, subid, subdir, sigdir, seg_name));

    %%% Parameters to convert graph to 3D skeleton
    %  sz =  the size of output volume. The order is [ y x z ]
    sz = size(seg);
    %  res = resolution of the nifti image [y,x,z] centimeters
    res = [0.0012, 0.0012, 0.0015];
    %  ds_flag = 1 or 0. (1=downsampled. 0=no-downsampled)
    ds_flag = 1;
    %  save_flag = to save the skeleton or not
    save_flag = 0;
    % Convert graph to skeleton
    [skel_pre] = sk3D(sz, Data.Graph, 'foo', res, ds_flag, save_flag);
    t_str = 'Unprocessed Volume';
    % Overlay the graph skeleton and the segmentation
    graph_seg_overlay(t_str, skel_pre, seg);
end

%% Plot graph and highlight nodes in loops
if visual
    % Copy edges into standard format
    s = edges(:,1); % source node
    t = edges(:,2); % target node
    % Create standard Matlab g
    g_mat = graph(s, t);
    % Find the nodes and edges belonging to loops
    [cnodes, cedges] = allcycles(g_mat);
    cnodes = unique(horzcat(cnodes{:}));
    cnodes = cnodes';
    
    %%% Create list of nodes to NOT down sample
    % Create array of node indices 
    nkeep = 1:1:size(nodes,1);
    nkeep = nkeep';
    % Find nodes not in loops
    nkeep = ~ismember(nkeep, cnodes);
    % Find array indices equal to 1. These are the node indices of non-loop
    % segments, which should be preserved during down sampling.
    nkeep = find(nkeep == 1);

    %%% Visualize graph
    visualize_graph(nodes, edges, 'Graph Before Loop Removal',nkeep) %#ok<*UNRCH> 
end

%% Remove Loops Parallel function (this calls the rm_loops for each loop)
% Move to mean minimum voxel intensity
v_min = 0.99;
% Down sample search radius
delta = 4;
% # iterations for mv2mean function in for-loop iteration in rm_loops
mv_iter = 1;
% Boolean for visualizing debugging graphs
viz = false;

rm_loop_tic = tic;
[~, ~] =...
    rm_loops_parallel_v0(nodes, edges, vol, delta, v_min, mv_iter, viz);
t_original = toc(rm_loop_tic);

rm_loop_tic = tic;
[~, ~] =...
    rm_loops_parallel(nodes, edges, vol, delta, v_min, mv_iter, viz);
t_improved = toc(rm_loop_tic);

tdiff = t_original - t_improved;
sprintf('The enhancement reduced overall computation time by %f second',...
    tdiff);