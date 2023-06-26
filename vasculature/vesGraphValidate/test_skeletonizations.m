%% Test the function downsample_nodes
%{
This script is for debugging the downsampling issue. The downsampling
currently connects disparate vessels. More details found here:
https://github.com/BUNPC/psoct_vessel_graphing/issues/1
%}
clear; clc; close all;
%% Initialize test bench
%%% Load vasculature binary mask
dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\AD_20832\dist_corrected\volume';
fname = 'ref_4ds_norm_inv_cropped_segment_sigma1_thresh0.25_crop.tif';
% Convert tif to mat
ves = TIFF2MAT(fullfile(dpath, fname));
% Convert to logical
ves = logical(ves);
% Voxel dimensions
vox_dim = [12, 12, 15];
% Minimum branch length for skeletonization and graphing
ves_min_len = 20;

%%% Convert mask to skeleton
ves_skel = bwskel(ves,'MinBranchLength', ves_min_len);
% Convert to skeleton to .TIFF for visualization
% tifout = 'ref_4ds_norm_inv_cropped_segment_sigma1_thresh0.25_skel.tif';
% ves_skel_im = im2uint8(ves_skel);
% segmat2tif(ves_skel_im, fullfile(dpath, tifout))

%% Convert skeleton to graph
% Minimum length of skeleton to convert to graph
skel_min_len = 20;

%%% Skel2Graph (new function)
% Minimum skeleton branch length for converting to graph
[a, node, link] = Skel2Graph3D(ves_skel, skel_min_len);

%%% seg_to_graph (old function)
[g2] = seg_to_graph(ves, vox_dim, ves_min_len);

%% Compare results
%%% Skel2Graph
% Convert to matlab graph
g1 = graph(a, 'omitselfloops');
% Copy node coordinates
x = vertcat(node.comx);
y = vertcat(node.comy);
z = vertcat(node.comz);
% Plot matlab graph struct
figure;
p1 = plot(g1,'XData',x,'YData',y,'ZData',z,'NodeColor','b');
title('Skel2Graph'); xlabel('x'); ylabel('y'); zlabel('z')
view(3)

%%% seg_to_graph
% Retrieve nodes/edges
nodes = g2.nodes;
edges = g2.edges;
% Copy edges into standard format
s = edges(:,1); % source node
t = edges(:,2); % target node
% Create standard Matlab graph
g2_mat = graph(s, t, 'omitselfloops');
% Plot matlab graph struct. This extracts the [x,y,z] of each node.
figure;
p2 = plot(g2_mat, 'XData', nodes(:,2), 'YData', nodes(:,1), 'ZData', nodes(:,3));
title('seg to graph'); xlabel('x'); ylabel('y'); zlabel('z')
view(3);
% Find cycles (loops)

%}

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





