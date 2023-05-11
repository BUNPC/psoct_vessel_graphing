%% This script is for testing the "rm_loops" function
% TODO: determine why the graph does not match the vessel segmentation.
%       the order of x,y,z in the nodes may be incorrect.

%% Graph struct from PSOCT graph.
clear; clc; close all;

% Load PSOCT graph
dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_human_brain_resources\test_data\Ann_Mckee_samples_10T\AD_20832\dist_corrected\volume';
fname = 'ref_4ds_norm_inv_cropped_segment_sigma1_thresh0.25_graph.mat';
load(fullfile(dpath, fname));
vox = Graph.vox;

% Downsample
% v = zeros(length(Graph.nodes),1);
% [nodes, edges,verifiedNodes,verifiedEdges] =...
%     downsample_nodes(Graph.nodes, Graph.edges, v, vox(1), vox(3));

nodes = Graph.nodes;
edges = Graph.edges;

% Copy edges into standard format
s = edges(:,1); % source node
t = edges(:,2); % target node

% Create standard Matlab graph
g = graph(s, t);

% Plot graph before removing loops
figure;
p = plot(g, 'XData', nodes(:,2), 'YData', nodes(:,1), 'ZData', nodes(:,3));
title('Graph Before Removing Loops'); xlabel('x'); ylabel('y'); zlabel('z')
view(3);

%% Find and delete cycles.
% cycles = node indices. edgecylces = edge indices
% [cycles,edgecycles] = allcycles(g, 'MaxNumCycles', 100, 'MaxCycleLength', 100);
[cycles,edgecycles] = allcycles(g);

% Highlight edges
for ii=1:length(edgecycles)
    highlight(p,'Edges',edgecycles{ii},'EdgeColor','r','LineWidth',1.5,'NodeColor','r','MarkerSize',6)
end

%%% Delete nodes that are part of cycles in graph
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
p = plot(H, 'XData', nodes(:,2), 'YData', nodes(:,1), 'ZData', nodes(:,3));
title('Graph after removing loops'); xlabel('x'); ylabel('y'); zlabel('z')
view(3)













