%% Test regraphNodes with list of nodes
%{
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

%% Initialize data path for linux or personal machine (debugging)
dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\';
% Subject IDs
subid = 'NC_6839';
subdir = '\dist_corrected\volume\';
sigdir = '\gsigma_1-3-5_gsize_5-13-21\';
% Segmentation filename
vol_name = 'ref_4ds_norm_inv_crop2.tif';
% Graph nodes/edges after marching ellipsoid
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
plot_graph(g_mat, nodes, 'Before Marching Ellipsoid')

%% Load list of nodes from multiple segments to regraph
fullpath = fullfile(dpath, subid, subdir, sigdir);
filename = strcat(fullpath, node_merge_struct);
g = load(filename);
regraph_idcs = g.node_merge_idx;
% Convert struct to cell
regraph_idcs = struct2cell(regraph_idcs);

%% Regraph nodes from group of segments
% Validated nodes
varray = zeros(length(nodes),1);
% Search distance (voxels)
delta = 4;
% Downsample
[segn, nodes, edges, validatedNodes,validatedEdges] =...
    regraphNodes_new(regraph_idcs, nodes, edges, varray, delta);

%% Call regraph for nodes from single segment

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