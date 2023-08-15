%% Test script for the CenterNodesXYZ function from David Boas
% David wrote code for moving nodes towards the center line of a segment.
% The original function was integrated into a GUI, but there is no
% documentation for the GUI. Therefore, this script will test a standalone
% version of the function.

% TODO: 
%   - wrap for-loop around mv_to_mean
%   - compare across iterations

clear; clc; close all;

%% Setup filename + processing variables
dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\';
subid = 'NC_6839';
subdir = '\dist_corrected\volume\';
sigdir = 'gsigma_1-3-5_gsize_5-13-21\';
vdata = 'ref_4ds_norm_inv_crop2.tif';
gdata = 'ref_4ds_norm_inv_crop2_segment_pmin_0.23_mask40_graph_data.mat';

% Imaging voxel intensity threshold
im_thresh = 52000;

%% Load PSOCT graph
tmp = load(fullfile(dpath, subid, subdir, sigdir, gdata), 'Data');
tmp = tmp.Data;

% Create new variable to match format of function
% im.angio = tmp.angio;
im.nodes = tmp.Graph.nodes;
im.edges = tmp.Graph.edges;

%% Load volumetric information
vol = TIFF2MAT(fullfile(dpath, subid, subdir, vdata));
im.angio = vol;

%% Test move to mean
im_mv = im;
for i=1:10
    im_mv = mv_to_mean(im_mv, im_thresh);
end

%%% Compare results
% Visualize uncentered graph
graph_vis(im.nodes, im.edges, 'Graph Before Moving to Mean')
% Visualize centered graph
graph_vis(im_mv.nodes, im_mv.edges, 'Graph After Moving to Mean')

%% Centering
% Unsure what this variable does
centerStep1vox = 0;
% Visualization
visualize_flag = 0;
% Run centering function
im_centered = center_nodes_xyz(im, centerStep1vox, visualize_flag, im_thresh);

%%% Compare results
% Visualize uncentered graph
graph_vis(im.nodes, im.edges, 'Graph Before Centering')
% Visualize centered graph
graph_vis(im_centered.nodes, im_centered.edges, 'Graph After Centering')

%% Visualize graph
function graph_vis(nodes, edges, title_str)
% Copy edges into standard format

s = edges(:,1); % source node
t = edges(:,2); % target node

% Create standard Matlab graph
g = graph(s, t);

% Plot graph before removing loops
figure;
p = plot(g, 'XData', nodes(:,1), 'YData', nodes(:,2), 'ZData', nodes(:,3));
title(title_str); xlabel('x'); ylabel('y'); zlabel('z')
view(3);
end