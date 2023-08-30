%% Main script for analyzing outpt of marching ellipsoid
%{
TODO:
- Revise regraph to only work along select segments
- Create a cylinder of radius r for the segment
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

%%% Visualization of graph
xlims = [0, 100]; ylims = [150, 300]; zlims = [50, 100];
xlims2 = [35,70]; ylims2 = [200, 260]; zlims2 = [70, 90];

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
plot_graph(g_mat, nodes, 'Before Marching Ellipsoid');
xlim(xlims2); ylim(ylims2); zlim(zlims2);

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
xlim(xlims2); ylim(ylims2); zlim(zlims2);

%% Regraph with old function (for comparison)
%%% Regraph
% Search delta
delta = 2;
% all nodes are unvalidated
validated_nodes = zeros(size(nodes,1),1);
% Regraph
[nodes, edges, ~,~] =...
    regraphNodes_new(nodes, edges,validated_nodes, delta);
% Convert to matlab graph struct
s = edges(:,1); % source node
t = edges(:,2); % target node
% Create graph
g_regraphed = graph(s, t);

%%% Plot graph with builtin function
tstr = {'Regraph (old function)','After Marching Ellipsoid'};
plot_graph(g_regraphed, nodes, tstr);
xlim(xlims2); ylim(ylims2); zlim(zlims2);

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
% Iterate over total number of segments
for ii = 1:max(bins)
    % Find number of nodes in segment (n)
    n = length(bins(bins==ii));
    % Convert number of nodes to node indices
    idcs = j : (j + n - 1);
    % Iterate counter
    j = max(idcs) + 1;
    % If <= 3 nodes, store n and node indices (this may not be used).
    if (1 <= n) && (n < nmin)
        del_nodes = [del_nodes, idcs];
    end

    % Store indices of nodes for segment
    nstruct(ii).node_idcs = idcs;
    
    % Find the index in the middle of the segment
    nmid = round(median(idcs));
    
    % Store the coordinates of the middle index. Later compare this to
    % the middle index of adjacent segments to determine if they are
    % within a distance threshold for combining.
    nstruct(ii).middle = nodes(nmid,:);
    
    % Find the normalized vector of the segment
    v = nodes(idcs(end),:) - nodes(idcs(1),:);
    nstruct(ii).v = v ./ norm(v);
    
    % TODO: Create a cylinder of radius r for the segment
end

%% Find segments w/ centers within dmin & norm(vectors) < threshold

%%% Identify nodes within dmin distance
% Convert struct to vertical array
x = vertcat(nstruct.middle);
% Create distance matrix (distance between each node)
d = squareform(pdist(x));
% Find index of node w/ euclidean distance below threshold
dmin = 30.0;
dmat = (d <= dmin);

%%% Identify segments w/ norm(vectors) < threshold
% Convert struct to vertical array
v = vertcat(nstruct.v);
% Calculate difference in vector of each node
vd = squareform(pdist(v));
% Find index of node w/ euclidean distance of vectors below threshold
vmin = 1.2;
vmat = (vd < vmin);

%%% Plot distrubution of difference matrices
figure; histogram(triu(d,1)); title('Coordinate Difference Matrix');
figure; histogram(triu(vd,1)); title('Vector Difference Matrix');

%%% Find segments satisfying both conditions
merge_mat = logical(dmat .* vmat');
merge_mat = triu(merge_mat, 1);
% Find nonzero elements in matrix
cc = bwconncomp(merge_mat,8);
% Find connected components
seg_merge_idx_struct = struct();
midx = 1;
for ii = 1 : cc.NumObjects
    % Convert from cell to array for ease of storage.
    cc_tmp = cc.PixelIdxList{ii};
    % If more than one pixel, then it's a group of connected pixels.
    % These are nodes that meet both conditions (distance and vector).
    % These will later be regraphed.
    if length(cc_tmp) > 1
        seg_merge_idx_struct(midx).idcs = cc_tmp;
        midx = midx + 1;
    end
end

%%% Convert node indices -> subscripts for groups meeting both conditions
% The indices within "seg_merge_idx_struct.idcs" reference the
% distance matrix. The row/column of the distance matrix are the node
% indices of nodes meeting both conditions. This section will:
%   - Convert the index to subscript (row, column)
%   - create an array of both the row and column subscripts
%   - Find the unique elements

group_node_idcs = struct();
% Convert segment indices to subscript
for ii = 1 : length(seg_merge_idx_struct)
    % Distance matrix indices of group of nodes meeting conditions
    merge_idx = seg_merge_idx_struct(ii).idcs;
    % Convert indices of segment to matrix subscripts
    [row, col] = ind2sub(size(vmat), merge_idx);
    % Create array with both row and column subscripts.
    % Only take the unique elements. 
    node_merge_idcs = unique(vertcat(row, col));
    % Extract nodes indices from all segments in group
    group_node_idcs(ii).node_idcs = horzcat(nstruct(node_merge_idcs).node_idcs);
end

%%% Save graph output
fout = '__adjacent_node_merge_list.mat';
fout = strcat(graph_name(1:end-4), fout);
node_merge_output = strcat(fullpath, fout);
save(node_merge_output, 'group_node_idcs');

%% Regraph the nodes within each segment
% Variables for downsampling
validatedNodes = zeros(size(nodes,1), 1);
delta = 2;
% Downsample (use new down sample function instead here)
[~, nodes_ds, edges_ds,~,~] = ...
    regraphNodes_new(group_node_idcs(ii).node_idcs,nodes,edges,validatedNodes,delta);
%%% Plot downsampled graph
s = edges_ds(:,1); % source node
t = edges_ds(:,2); % target node
g_me_ds = graph(s, t);
plot_graph(g_me_ds, nodes_ds, 'Marching Ellipsoid - Regraphed')

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
% Downsample (use new down sample function instead here)
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

%% Function to plot graph from matlab graph
function plot_graph(g, nodes, title_str)
figure;
p = plot(g, 'XData', nodes(:,1), 'YData', nodes(:,2), 'ZData', nodes(:,3));
p.NodeColor = 'red';
p.EdgeColor = 'blue'; p.LineWidth = 2;

xlabel('x'); ylabel('y'); zlabel('z'); title(title_str);
set(gca, 'FontSize', 25); grid on;
view(3);
end

%% Plot graph with lines
function scatter_graph(edges, nodes)
% Scatter plot from edges and nodes
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
end
