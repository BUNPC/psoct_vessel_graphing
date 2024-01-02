%% Purpose: overlay segmentation, skeleton, and graph
%{
Purpose:
- load a graph structure
- Call remove loops function
%}
% clear; close all; clc;

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

%% Initialize data paths
if ispc
    % Top-level directories
    dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\';
    subid = 'AD_21424';
%     subid = 'NC_21499';
%     subid = 'NC_6974';
    subdir = '\dist_corrected\volume\combined_segs';
    sigdir = 'gsigma_3-5-7_5-7-9_7-9-11\';
    % Segment and graph data
    seg_name = 'seg_masked.mat';
    graph_name = 'seg_maskedgraph_data.mat';
    
    %%% Output Data filenames
    skel_out = strcat(graph_name(1:end-4),'_loops_rm.tif');
    skel_out = char(fullfile(dpath, subid, subdir, sigdir, skel_out));
    gdata_out = strcat(graph_name(1:end-4),'_loops_rm.mat');
    
elseif isunix
    %%% CTE_6489
    % Top-level directories
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
    subid = 'CTE_6912';
    subdir = '/dist_corrected/volume/';
    sigdir = 'gsigma_1-2-3-4-5_gsize_5--9-13-17-21/';
    vdata = 'ref_4ds_norm_inv.tif';
    % Data with many nested loops (probability threshold = 0.21)
    seg_name = 'ref_4ds_norm_inv_segment_pmin_0.26_mask40.tif';
    graph_name = 'ref_4ds_norm_inv_segment_pmin_0.26_mask_40_graph_data.mat';
    %}
    
    %%% Output Data filenames
    skel_out = strcat(graph_name(1:end-4),'_loops_rm.tif');
    skel_out = char(fullfile(dpath, subid, subdir, sigdir, skel_out));
    gdata_out = strcat(graph_name(1:end-4),'_loops_rm.mat');
end

%% Import graph and segmentation
% Define entire filepath 
fullpath = fullfile(dpath, subid, subdir, sigdir);
% Import segmentation
fname = fullfile(fullpath, seg_name);
segm = load(fname, 'segm');
segm = segm.segm;
% Import graph data
fname = fullfile(fullpath, graph_name);
data = load(fname, 'Data');
data = data.Data;
graph = data.Graph;

%% Convert graph -> skel. Overlay skeleton and segmentation.
%%% Determine size of output file
% Size of matrix with segmentation [x,y,z]
sz = size(segm);

%%% Parameters to convert graph to 3D skeleton
% Dummy variable for resolution of NIFTI. This won't be used, since the
% save nifti flag will be false.
res = [0.01, 0.01, 0.01];
% Down sample flag. Set true if the graph was down sampled
ds_flag = 0;
% Figure title
fig_name = subid;

%%% Overlay skel + graph
[seg, skel] = overlay_seg_skel_graph(sz, graph, res, ds_flag, fig_name, segm);

%% Overlay graph and segmentation

%%% Convert segmentation to graph (nodes and edges)
g = seg_to_graph(seg, [12,12,15]);
visualize_graph(g.nodes, g.edges, 'Subset',[])

%%% Call regraph (downsample)
validated_nodes = zeros(1,length(g.nodes));
delta = 2;
[g.nodes, g.edges, ~, ~] =...
    regraphNodes_new(g.nodes, g.edges, validated_nodes, delta);

%%% Call Move to Mean
% Add segmentation to struct
seg = permute(seg, [2,1,3]);
g.angio = seg;
% Call move to mean
g2 = mv_to_mean(g, 0.1);

%% Visualize graph
scatter_graph(g2.edges, g2.nodes, 'Downsampled Subset'); grid off;

%% Save the graph into the format for vesGraphValidat
Data = struct();
Data.angio = seg;
Data.Graph.edges = g.edges;
Data.Graph.nodes = g.nodes;
fout = fullfile(fullpath, 'cropped_graph.mat');
save(fout,'Data','-v7.3');








