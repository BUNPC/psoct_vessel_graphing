%% Test the function downsample_nodes
%{
This script is for debugging the downsampling issue. The downsampling
currently connects disparate vessels. More details found here:
https://github.com/BUNPC/psoct_vessel_graphing/issues/1
%}
clear; clc; close all;
%% Initialize test bench
%%% Load test data set
dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\AD_20832\dist_corrected\volume';
fname = 'ref_4ds_norm_inv_cropped_segment_sigma1_thresh0.25_graph.mat';
load(fullfile(dpath, fname));
% Voxel resolution [x, y, z] (microns)
vox = Graph.vox;
% Retrieve nodes/edges
nodes = Graph.nodes;
edges = Graph.edges;

%%% Create matlab graph data structure.
% Copy edges into standard format
s = edges(:,1); % source node
t = edges(:,2); % target node
% Create standard Matlab graph
g = graph(s, t);
% Plot graph. This extracts the [x,y,z] of each node.
figure;
p = plot(g, 'XData', nodes(:,2), 'YData', nodes(:,1), 'ZData', nodes(:,3));
title('Graph Before Downsample'); xlabel('x'); ylabel('y'); zlabel('z')
view(3);

%% Downsample w/ new matlab function
v = zeros(length(Graph.nodes),1);
[nodes_ds, edges_ds, ~, ~] =...
    downsample_graph(nodes, edges, vox(1), vox(3));

%%% Create matlab graph data structure.
% Copy edges into standard format
s = edges_ds(:,1); % source node
t = edges_ds(:,2); % target node
% Create standard Matlab graph
g_ds = graph(s, t);
% Plot graph. This extracts the [x,y,z] of each node.
figure;
p_ds = plot(g_ds, 'XData', nodes_ds(:,2), 'YData', nodes_ds(:,1), 'ZData', nodes_ds(:,3));
title('Graph After Downsample'); xlabel('x'); ylabel('y'); zlabel('z')
view(3);

% Ensure # segments is equal before/after downsampling

% Retrieve new nodes/edges
% nodes_ds = p.





