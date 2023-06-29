%% Debug the "rm_loops" function
%{
TODO:
- test with graph from entire volume + determine run time

%}


%% Create graph from PSOCT mask.
clear; clc; close all;

% Load PSOCT graph
dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\AD_20832\dist_corrected\volume\';
fname = 'ref_4ds_norm_inv_cropped_segment_sigma1_thresh0.25_crop2_pruned.tif';
ves_mask = TIFF2MAT(fullfile(dpath, fname));
vox_dim = [12, 12, 15];
min_branch_len = 20;

% Convert mask to graph
g = seg_to_graph(ves_mask, vox_dim);

% Retrieve nodes/edges
vox = g.vox;
nodes = g.nodes;
edges = g.edges;

% Copy edges into standard format
s = edges(:,1); % source node
t = edges(:,2); % target node

% Create standard Matlab graph
g = graph(s, t);

% Plot graph before removing loops
figure;
p = plot(g, 'XData', nodes(:,1), 'YData', nodes(:,2), 'ZData', nodes(:,3));
title('Graph Before Removing Loops'); xlabel('x'); ylabel('y'); zlabel('z')
view(3);

%% Find and delete cycles.
% cycles = node indices. edgecylces = edge indices
% [cycles,edgecycles] = allcycles(g, 'MaxNumCycles', 100, 'MaxCycleLength', 100);
[cycles,edgecycles] = allcycles(g);

%%% Highlight edge cycles
for ii=1:length(edgecycles)
    highlight(p,'Edges',edgecycles{ii},'EdgeColor','r','LineWidth',1.5,'NodeColor','r','MarkerSize',6)
end

%%% Highlight cycles
for ii=1:length(cycles)
    highlight(p,'Edges',cycles{ii},'EdgeColor','g','LineWidth',1.5,'NodeColor','g','MarkerSize',6)
end

%%% Matrix of cycles in graph
% Convert cell to matrix
cyclical = [];
for i=1:length(cycles)
    % Add to matrix
    cyclical = [cyclical, cycles{i}];
end
% Find unique nodes
u = unique(cyclical);

%%% Remove cyclical nodes
H = rmnode(g, u);
% Remove corresponding entries in x,y,z coordinates
nodes(u,:) = [];

%%% Plot updated graph
figure;
p = plot(H, 'XData', nodes(:,1), 'YData', nodes(:,2), 'ZData', nodes(:,3));
title('Graph after removing loops'); xlabel('x'); ylabel('y'); zlabel('z')
view(3)

%% Find intersection of cycle and segment
%{
Deleting the cycle results in a discontinuity in the vessel. To ensure a
continuous vessel, we must first locate the intersections between the cycle
and the segment, then we must create an edge between these points.
%}

% Find connected components of binary vessel mask
ves_skl = bwskel(logical(ves_mask),'MinBranchLength', min_branch_len);   
cc = bwconncomp(ves_skl);

%%% Find branch points within connected component (vessel)
br = bwmorph3(ves_skl, 'branchpoints');
% Find matrix indices of branch points
brix = find(br);

%%% Find 2 branch points (n1, n2) belonging to the same cycle
% Struct for storing branch points for each segment
bruct = struct();

%%% Iterate over cycles (or groups, or edgecycles?)
for ii = 1:length(cycles)
    % Find branchpoint indices within the group
    ctmp = cycles{ii};
    idx = find(ismember(brix, ctmp));
    % Convert from brix index to ves_skl matrix index
    bruct(ii).idx = brix(idx);
end

% Remove cycles from graph

% Connect the two branch points that are in the same vessel












