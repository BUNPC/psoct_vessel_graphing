%% Main script for performing marching ellipsoid
%{
This package was initially created by collaborators of David Boas. It has
been modified over the years. This main script is still a work in progress.

This script was run with a subset of a larger g. The code did not
compute the graph2.edges. Need to debug this.
%}
clear; clc; close all;
%% Add top-level directory of code repository to path
% This allows Matlab to find the functions in the project folders

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

%%% Local machine
dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\';
% Subject IDs
subid = 'NC_6839';
subdir = '\dist_corrected\volume\';
sigdir = '\gsigma_1-3-5_gsize_5-13-21\';
% Segmentation filename
vol_name = 'ref_4ds_norm_inv_crop2.tif';
% Graph filename
graph_name = 'ref_4ds_norm_inv_crop2_segment_pmin_0.23_mask40_ds_mean_ds_graph.mat';
graph_me_name = 'ref_4ds_norm_inv_crop2_segment_pmin_0.23_mask40_ds_mean_ds_marching_ellipse.mat';
% graph_me_name = 'ref_4ds_norm_inv_crop2_segment_pmin_0.23_mask40_ds_mean_ds_mstep4_cutoff4';
% View flag for marching ellipsoid (1 = view).
viewflag = 1;

%% Load volume and g from Data
vox_dim = [12, 12, 15];

%%% Load volume
fullpath = fullfile(dpath, subid, subdir);
filename = strcat(fullpath, vol_name);
vol = mat2gray(TIFF2MAT(filename));
% volumeViewer(vol);

%%% Load graph
fullpath = fullfile(dpath, subid, subdir, sigdir);
filename = strcat(fullpath, graph_name);
g = load(filename);
g = g.im_re;
nodes = g.nodes;
edges = g.edges;
% Copy edges into standard format
s = edges(:,1); % source node
t = edges(:,2); % target node

% Create standard Matlab g
g_mat = graph(s, t);
plot_graph(g_mat, nodes, 'Before Marching Ellipsoid')

%% Load the graph after marching ellipsoid
fullpath = fullfile(dpath, subid, subdir, sigdir);
filename = strcat(fullpath, graph_me_name);
g = load(filename);
g = g.g;
nodes = g.nodes;
% Debugging line
edges = g.edges;
% Copy edges into standard format
s = edges(:,1); % source node
t = edges(:,2); % target node
% Create graph
g_mat = graph(s, t);

%%% Plot graph with builtin function
tstr = 'After Marching Ellipsoid';
plot_graph(g_mat, nodes, tstr);
% xlim([60, 140]); ylim([110, 180]); zlim([0,90]);

%%% Plot graph with lines
%{
figure;
% Scatter plot of nodes
scatter3(nodes(:,1), nodes(:,2), nodes(:,3), '.','b');
% Plot edges
for ii=1:length(edges)
    n1 = edges(ii,1);
    n2 = edges(ii,2);
    x = [nodes(n1,1), nodes(n2,1)];
    y = [nodes(n1,2), nodes(n2,2)];
    z = [nodes(n1,3), nodes(n2,3)];
    line(x, y, z, 'Color', 'red');
end
%}

%% Smooth results after marching ellipsoid
% Variables for downsampling
validatedNodes = zeros(size(nodes,1), 1);
delta = 4;
% Downsample
[~, nodes_ds, edges_ds,~,~] = ...
    regraphNodes_new([],nodes,edges,validatedNodes,delta);
%%% Plot downsampled graph
s = edges_ds(:,1); % source node
t = edges_ds(:,2); % target node
g_me_ds = graph(s, t);
plot_graph(g_me_ds, nodes_ds, 'Marching Ellipsoid - Regraphed')

%% Find segments with fewer than nmin nodes and find node indices
% Connectivity analysis: find which segment each node belongs to. This will
% output an array of indices (bins). Each entry in bins corresponds to the
% node index. The value of the entry corresponds to the segment ID. For
% example, if bins(1:3)==1, then nodes 1-3 belong to segment 1.
bins = conncomp(g_mat);
% Array to track nodes for deletion
del_nodes = [];
% Minimum number of nodes per edge to keep
nmin = 4;
% Track node index
j = 1;
% Struct to store node indices belonging to edges w/ fewer than nmin nodes
nstruct = struct();
% Index for storing data in nstruct
seg_idx = 1;
% variable to track number of nodes/segment
nmax = 1;
% Iterate over total number of segments
for ii = 1:max(bins)
    % Find number of nodes in segment (n)
    n = length(bins(bins==ii));
    % Convert number of nodes to node indices
    idcs = j : (j + n - 1);
    % Iterate counter
    j = max(idcs) + 1;
    % If <= 3 nodes, store n and node indices
    if (1 <= n) && (n < nmin)
        del_nodes = [del_nodes, idcs];
    else
        % Store indices of nodes
        nstruct(seg_idx).node_idcs = idcs;
        
        % Find the index in the middle of the segment
        nmid = round(median(idcs));
        
        % Store the coordinates of the middle index. Later compare this to
        % the middle index of adjacent segments to determine if they are
        % within a distance threshold for combining.
        nstruct(seg_idx).middle = nodes(nmid,:);
        
        % Find the normalized vector of the segment
        v = nodes(idcs(end),:) - nodes(idcs(1),:);
        nstruct(seg_idx).v = v ./ norm(v);
        
        % Create a cylinder of radius r for the segment
        
        % Iterate the index for the struct
        seg_idx = seg_idx + 1;
    end
end

%% Delete nodes/edges belonging to segments with < nmin nodes
%{
% Delete nodes from graph struct
g_mat_rm = rmnode(g_mat, del_nodes);
% Delete from the variables "nodes" and "edges"
nodes_rm = nodes;
nodes_rm(del_nodes,:) = [];
% Update edges
edges_rm = g_mat_rm.Edges;
edges_rm = edges_rm{:,:};
% Plot graph (minus segments w/ <= 3 nodes)
tstr = strcat('Minimum nodes/edge = ', num2str(nmin));
plot_graph(g_mat_rm, nodes_rm, {'Marching Ellipsoid', tstr})

%%% Downsample after removing segments
% Variables for downsampling
validatedNodes = zeros(size(nodes,1), 1);
delta = 8;
% Downsample
[~, nodes_ds, edges_ds,~,~] = ...
    regraphNodes_new([],nodes_rm,edges_rm,validatedNodes,delta);

%%% Create graph
s = edges_ds(:,1); % source node
t = edges_ds(:,2); % target node
g_me_ds = graph(s, t);

%%% Plot downsampled graph
plot_graph(g_me_ds, nodes_ds, {'After Marching Ellipsoid','Removing small segments & Regraphed'})
% xlim([60, 140]); ylim([110, 180]); zlim([0,90]);
%}

%% Merge segments w/ centers within dmin & norm(vectors) < threshold

%%% Identify nodes within dmin distance
% Create distance matrix
x = vertcat(nstruct.middle);
z = squareform(pdist(x));
% Find index of node w/ euclidean distance below threshold
dmin = 30.0;
dmat = (z <= dmin);
% didx = find( (z ~= 0) & (z < dmin));

%%% Identify segments w/ norm(vectors) < threshold
% Create distance matrix
x = vertcat(nstruct.v);
z = squareform(pdist(x));
% Find index of node w/ euclidean distance of vectors below threshold
vmin = 1.0;
vmat = (z < vmin);
vidx = unique(find( (z ~= 0) & (z < vmin)));

%%% Find nodes satisfying both conditions
merge_mat = dmat .* vmat';
merge_mat = triu(merge_mat, 1);

%%% TODO:
% find groups of connected components
% regraph edges within each group

%%% Find indices of nodes that satisfy both conditions
% Unique indices
merge_idx = unique(intersect(didx, vidx));
% Convert indices to matrix subscripts to find similar segments
[row, col] = ind2sub(size(z), merge_idx);

%%% Remove one of the redundant segments meeting conditions
% Create new node + edge list
node_merge = nodes;
edge_merge = edges;
% Initialize array for storing indices to delete
ndel = [];
for ii = 1:length(col)
    % Remove one of the duplicate segments (either row or column index)
    idx = col(ii);
    % Add indices to delete
    ndel = [ndel, nstruct(idx).node_idcs];
end
% Remove redundant nodes
nodes(ndel,:) = [];
% Remove edges containing node indices in ndel

%}

%% Function to plot graph from matlab graph
function plot_graph(g, nodes, title_str)
figure;
p = plot(g, 'XData', nodes(:,1), 'YData', nodes(:,2), 'ZData', nodes(:,3));
p.EdgeColor = 'red'; p.LineWidth = 1.5;
xlabel('x'); ylabel('y'); zlabel('z'); title(title_str);
set(gca, 'FontSize', 20); grid on;
view(3);
end