%% Test regraphNodes with list of nodes
%{
Purpose:
- load a graph structure
- Call remove loops function
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
if ispc
    %%% AD 20832
    % Top-level directories
    dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\';
    subid = 'AD_20832';
    subdir = '\dist_corrected\volume\';
    sigdir = 'gsigma_1-3-5_gsize_5-13-21\';
    vdata = 'ref_4ds_norm_inv_crop_small.tif';
    % Segment and graph data
    seg_name = 'ref_4ds_norm_inv_crop_small_segment_pmin_0.26.tif';
    gdata = 'ref_4ds_norm_inv_crop_small_segment_pmin_0.26_graph_data.mat';
    
    %%% NC_6839
    %{
    % Top-level directories
    dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\';
    subid = 'NC_6839';
    subdir = '\dist_corrected\volume\';
    sigdir = 'gsigma_1-3-5_gsize_5-13-21\';
    vdata = 'ref_4ds_norm_inv_crop2.tif';
    % Data with many nested loops (probability threshold = 0.21)
    seg_name = 'ref_4ds_norm_inv_crop2_segment_pmin_0.21.tif';
    gdata = 'ref_4ds_norm_inv_crop2_segment_pmin_0.21_mask_40_graph_data.mat';
    %}

    %%% Test Data from Etienne
    %{
    % Top-level directories
    dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\Martinos_Datasets\data_stats\version_8\post_processed';
    %}

    %%% Output Data filenames
    skel_out = strcat(gdata(1:end-4),'_loops_rm.tif');
    skel_out = char(fullfile(dpath, subid, subdir, sigdir, skel_out));
    gdata_out = strcat(gdata(1:end-4),'_loops_rm.mat');
    
elseif isunix
    %%% AD 20832
    %{
    % Top-level directories
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
    subid = 'AD_20832';
    subdir = '/dist_corrected/volume/';
    sigdir = 'gsigma_1-3-5_gsize_5-13-21/';
    vdata = 'ref_4ds_norm_inv_crop_small.tif';
    % Segment and graph data
    seg_name = 'ref_4ds_norm_inv_crop_small_segment_pmin_0.26.tif';
    gdata = 'ref_4ds_norm_inv_crop_small_segment_pmin_0.26_graph_data.mat';
    %}
    
    %%% NC_6839
    %{
    % Top-level directories
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
    subid = 'NC_6839';
    subdir = '/dist_corrected/volume/';
    sigdir = 'gsigma_1-3-5_gsize_5-13-21/';
    vdata = 'ref_4ds_norm_inv_crop2.tif';
    % Data with many nested loops (probability threshold = 0.21)
    seg_name = 'ref_4ds_norm_inv_crop2_segment_pmin_0.21.tif';
    gdata = 'ref_4ds_norm_inv_crop2_segment_pmin_0.21_mask_40_graph_data.mat';
    %}

    %%% CTE_6489
    % Top-level directories
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
    subid = 'CTE_6912';
    subdir = '/dist_corrected/volume/';
    sigdir = 'gsigma_1-2-3-4-5_gsize_5--9-13-17-21/';
    vdata = 'ref_4ds_norm_inv.tif';
    % Data with many nested loops (probability threshold = 0.21)
    seg_name = 'ref_4ds_norm_inv_segment_pmin_0.26_mask40.tif';
    gdata = 'ref_4ds_norm_inv_segment_pmin_0.26_mask_40_graph_data.mat';
    %}
    

    %%% Test Data from Etienne
    %{
    % Top-level directories
    dpath = '/projectnb/npbssmic/ns/martinos_etienne_test_data/data_stats/version_8/post_processed/';
    % TODO:
    %   - read in the NIFTI file
    %   - verify same skel_out filepaths
    mri = MRIread(strcat(dpath, fname));
    fprintf('mri loaded.\n')
    mri.vol = single(mri.vol);
    %}

    %%% Output Data filenames
    skel_out = strcat(gdata(1:end-4),'_loops_rm.tif');
    skel_out = char(fullfile(dpath, subid, subdir, sigdir, skel_out));
    gdata_out = strcat(gdata(1:end-4),'_loops_rm.mat');
end
%% Load PSOCT graph, volume, segmentation

%%% Load Graph
Data = load(fullfile(dpath, subid, subdir, sigdir, gdata), 'Data');
Data = Data.Data;
nodes = Data.Graph.nodes;
edges = Data.Graph.edges;

%%% Load volumetric information and set threshold
% Import volume
vol = TIFF2MAT(fullfile(dpath, subid, subdir, vdata));
% Imaging voxel intensity threshold (normalized [0,1])
im_thresh = 0.75;
% Check for data type
if strcmp(class(vol),'uint16')
    r = 65535;
    im_thresh = im_thresh .* r;
elseif strcmp(class(vol),'uint8')
    r = 255;
    im_thresh = im_thresh .* r;
end

%%% Load segmentation stack
% Import volume
seg = TIFF2MAT(fullfile(dpath, subid, subdir, sigdir, seg_name));
if visual == true
    volshow(seg);
end

%%% Parameters to convert graph to 3D skeleton
%  sz =  the size of output volume. The order is [ y x z ]
sz = size(seg);
%  res = resolution of the nifti image [y,x,z] centimeters
res = [0.0012, 0.0012, 0.0015];
%  ds_flag = 1 or 0. (1=downsampled. 0=no-downsampled)
ds_flag = 1;
%  save_flag = to save the skeleton or not
save_flag = 0;

%%% Overlay graph and segmentation
% Convert graph to skeleton
[skel_pre] = sk3D(sz, Data.Graph, 'foo', res, ds_flag, save_flag);
t_str = 'Unprocessed Volume';
% Overlay the graph skeleton and the segmentation
if visual
    graph_seg_overlay(t_str, skel_pre, seg);
end

%%% Create list of nodes/edges in loops to down sample
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
if visual
    visualize_graph(nodes, edges, 'Graph Before Loop Removal',nkeep)
end
%% Call Remove Loops function (move2mean and down sample loops)

% Move to mean minimum voxel intensity
v_min = 0.99;

% Down sample search radius
delta = 6;

% # iterations for mv2mean function in for-loop iteration in rm_loops
mv_iter = 1;

% boolean for whether or not to visualize debugging graph figures
viz = true;

[node_rm, edges_rm] =...
    rm_loops(nodes, edges, vol, delta, v_min, mv_iter, viz);

%% Visualize Results and Save

%%% Overlay the result with the segmentation
% Assign nodes/edges to graph for creating skeleton
g.nodes = node_rm;
g.edges = edges_rm;
% Create skeleton of graph
[skel] = sk3D(sz, g, 'foo', res, ds_flag, save_flag);
% Save output skeleton
segmat2tif(skel, skel_out);

% Overlay skeleton with segmentation
if visual
    % Display graph with all nodes green
    nidx = 1:size(node_rm,1);
    visualize_graph(node_rm, edges_rm, 'Graph After Loop Removal',nidx);
    % Overlay skeleton with segmentation
    fig_title = 'Loops Removed';
    graph_seg_overlay(fig_title, skel, seg)
end

%%% Save outputed graph
save(fullfile(dpath, subid, subdir, sigdir, gdata_out),'g','-v7.3');

%% Overlay graph and segmentation
function graph_seg_overlay(fig_title, skel, seg)
%graph_seg_overlay Overlay the graph and segmentation to visually verify
% INPUTS:
%   fig_title (string): name of figure
%   skel (matrix): skeleton of graph (binary)
%   seg (matrix): segmentation (binary)

% Initialize the 3D figure properties
view_panel = uifigure(figure,'Name',fig_title); close;
v = viewer3d(view_panel);
v.BackgroundColor = 'w';
v.BackgroundGradient = 'off';

% Display volume of skeleton (from graph)
h = volshow(skel,'Parent',v);
h.Parent.BackgroundColor = 'w';

% Overlay the segmentation
h.OverlayData = seg;
h.OverlayAlphamap = 0.1;
end