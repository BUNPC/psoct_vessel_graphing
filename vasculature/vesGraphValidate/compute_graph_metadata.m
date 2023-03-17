%% Generate metadata for Graph.[nodes, edges]
% This script calls the function nodeGrps_vesSegment.m, which determines
% the metadata below:
% im.nodeGrp = nodeGrp;
% im.nodeSegN = nodeSegN;
% im.segNedges = segNedges;
% im.segLen = segLen;
% im.segLen_um = segLen_um;
% im.edgeSegN = edgeSegN;
% im.segEndNodes = segEndNodes;
% im.segPos = squeeze(mean(reshape(nodePos(im.segEndNodes,:),[2 length(segLen) 3]),1));

%% Setup Environment
% Top-level directory
addpath '\\ad\eng\users\m\h\mhyman\My Documents\Boas_Lab\psoct_human_brain';
% Skeleton recon directory
addpath '\\ad\eng\users\m\h\mhyman\My Documents\Boas_Lab\psoct_human_brain\vasculature\Vessel_reconstruction_graph_refinement\Skeleton';

% Hardcoded file for debugging:
pathname = fullfile(pwd, '..\');    % contained in subfolder
filename = 'frangi_seg.mat';

%% Load the Graph struct
tmp = load([pathname, filename]);
nodes = tmp.Graph.nodes;
edges = tmp.Graph.edges;

%% Call the metadata function
md = nodeGrps_vesSegment(nodes, edges);

%% copy metadata to the Graph struct
for fn = fieldnames(md)'
   Graph.segInfo.(fn{1}) = md.(fn{1});
end

%% Add angio matrix to Data
% clear; 
% load('Data.mat')

% Load the I_seg matrix from the vesSegment.m function
pathname = fullfile(pwd, '..\vesSegment\');    % contained in subfolder
filename = 'ves_seg.mat';
angio = load([pathname, filename], 'I_seg');
angio = angio.I_seg;

%% Save output
save('Data2', 'Graph', 'angio');