%% Test regraphNodes with list of nodes
%{
This test script was developed for downsampling a subset of nodes. The goal
was to connect parallel segments but not antiparallel ones.

The purpose of this test script is to test the following functionality of
the regraphNodes_new function:
- iterate over a list of nodes
- list of nodes from one segment
- list of nodes from multiple segments

Overview:
- load a graph structure
- call regraph for nodes from single segment
- call regraph for nodes from adjacent segments
%}

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
topdir = mydir(1:idcs(end));
addpath(genpath(topdir));

%% Initialize data paths
dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\';
% Subject IDs
subid = 'NC_6839';
subdir = '\dist_corrected\volume\';
sigdir = '\gsigma_1-3-5_gsize_5-13-21\';
% PS-OCT volume filename
vol_name = 'ref_4ds_norm_inv_crop2.tif';
% Struct of nodes/edges after applying the marching ellipsoid
graph_me_name = 'ref_4ds_norm_inv_crop2_segment_pmin_0.23_mask40_ds_mean_ds_marching_ellipse.mat';
% List of node indices within each list of segments to merge
node_merge_struct = 'ref_4ds_norm_inv_crop2_segment_pmin_0.23_mask40_ds_mean_ds_graph__adjacent_node_merge_list';

%% Load graph nodes/edges from marching ellipsoid graph

fullpath = fullfile(dpath, subid, subdir, sigdir);
filename = strcat(fullpath, graph_me_name);
g = load(filename);
g = g.g;
nodes = g.nodes;
edges = g.edges;
% Copy edges into standard format
s = edges(:,1); % source node
t = edges(:,2); % target node
% Create standard Matlab g
g_mat = graph(s, t);
% plot_graph(g_mat, nodes, 'After Marching Ellipsoid (unfiltered)')

%% Load list of nodes to down sample
% This struct contains a field (node_merge_idx). The rows in this field
% correspond to groups of nodes meeting the downsampling criteria from the
% marching ellipsoid output. These node indices will be downsampled.

% Create filepath
fullpath = fullfile(dpath, subid, subdir, sigdir);
filename = strcat(fullpath, node_merge_struct);
% Load struct
node_ds_idcs = load(filename);
node_ds_idcs = node_ds_idcs.group_node_idcs;
% Convert struct to cell
node_ds_idcs = struct2cell(node_ds_idcs);

%% Down sample nodes from group of segments
% Search distance (voxels)
delta = 4;

% Downsample with old method for debugging
validatedNodes = zeros(size(nodes,1),1);
[nodes_old_method, edges_old_method, ~,~] = ...
    regraphNodes_new(nodes, edges, validatedNodes, delta, delta);

%% Downsample with new method
[nodes_ds, edges_ds] =...
    downsample_subset_me(node_ds_idcs, nodes, edges, delta);



%%% Plot with scatterplot and lines
graph_title_str = 'Down Sampled Marching Ellipsoid Segment Groups';
scatter_graph(edges_ds, nodes_ds, graph_title_str);

%%% Plot result
% Copy edges into standard format
s = edges_ds(:,1); % source node
t = edges_ds(:,2); % target node
% Create standard Matlab g
g = graph(s, t);
graph_title_str = 'Down Sampled Marching Ellipsoid Segment Groups';
% Plot matlab graph
plot_graph(g, nodes_ds, graph_title_str);

%% Call regraph for nodes from entire graph
% Validated nodes
varray = zeros(length(nodes),1);
[nodes, edges, validatedNodes,validatedEdges] =...
    regraphNodes_new(nodes, edges, varray, delta, delta);

%% Create list of nodes for each segment

%% Function to plot graph from matlab graph
function plot_graph(g, nodes, title_str)
figure;
p = plot(g, 'XData', nodes(:,1), 'YData', nodes(:,2), 'ZData', nodes(:,3));
p.EdgeColor = 'red'; p.LineWidth = 1.5;
xlabel('x'); ylabel('y'); zlabel('z'); title(title_str);
set(gca, 'FontSize', 20); grid on;
view(3);
end