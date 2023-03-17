%% Remove the spur segments
%{
This script uses a nested struct "Data.Graph" which contains metadata about
the Graph struct. The Graph struct is created in the script SegtoGraph.m.

TODO:
- Determine meaning of A_idx_end{1,1}
- Determine logic of last two steps (new edges and new nodes)
%}

%% Prep environment & load Data, 
% clear; close all; clc;
% tmp = load('Data_processed.mat');
% Data = tmp.Output;
global Data

%% Parse segments
% Assign local variables for ease of access
nodes = Data.Graph.nodes;
edges = Data.Graph.edges;
% Find the segments <= 10 units in length
seg_del = find(Data.Graph.segInfo.segLen <= 10);
% Find intersections between [segments <= 10] & A_idx_end{1,1}
seg_del = intersect(Data.Graph.segInfo.A_idx_end{1,1}, seg_del);

% Iterate through segments and find indices
idx_del = [];
for i=1:length(seg_del)
    idx =  find(Data.Graph.segInfo.nodeSegN == seg_del(i));
    endnodes = Data.Graph.segInfo.segEndNodes(seg_del(i),:);
    endnodes = endnodes(:);
    segs1 = find(Data.Graph.segInfo.segEndNodes(:,1) == endnodes(1) | Data.Graph.segInfo.segEndNodes(:,2) == endnodes(1));
    segs2 = find(Data.Graph.segInfo.segEndNodes(:,1) == endnodes(2) | Data.Graph.segInfo.segEndNodes(:,2) == endnodes(2));
    
    % TODO: both of the if statements may be run, in which case idx would
    % be set twice. Do we need both if statements?
    if length(segs1) > 1
        tsegs1 = setdiff(segs1,seg_del);
        if ~isempty(tsegs1)
            idx = setdiff(idx,endnodes(1));
        end
    end
    if length(segs2) > 1
        tsegs2 = setdiff(segs2,seg_del);
        if ~isempty(tsegs2)
            idx = setdiff(idx,endnodes(2));
        end
    end
    % Add index to array
    idx_del = [idx_del; idx];
end

%% Create nodeMap, which is used to create an updated set of edges (purpose?)
% variable "nodeMap" should be renamed to "edgeMap"

nNodes = size(nodes,1);      % nNodes = total number of nodes
map = (1:nNodes)';           % map = [1, 2, ..., nNodes]
map(idx_del) = [];           % Remove array entries corresponding to idx_del
mapTemp = (1:length(map))';  % mapTemp = [1, 2, ..., length(map)]
nodeMap = zeros(nNodes,1);   % nodeMap = [0, 0,...,0]
nodeMap(map) = mapTemp;      % assign nonzero elements from map --> nodeMap

% Remove zero values from edges array.
edgesNew = nodeMap(edges);
zero_idx = find(edgesNew(:,1) == 0 | edgesNew(:,2)==0);
edgesNew(zero_idx,:) = [];

%% Remove indices from nodes
nodesNew = nodes;
nodesNew(idx_del,:) = [];

%% Save graph
Data.Graph.nodes = nodesNew;
Data.Graph.edges = edgesNew;
Data.Graph.verifiedNodes = zeros(size(Data.Graph.nodes,1),1);
disp([num2str(length(idx_del)) ' segments detected and removed from graph']);

save('Data_segpruned.mat', 'Data');