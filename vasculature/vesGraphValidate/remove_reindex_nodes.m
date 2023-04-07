function [pos_reindex, edge_reindex] =...
    remove_reindex_nodes(rm_list, nodePos, nodeEdges)
% Remove nodes from nodePos and nodeEdges. Reindex nodeEdges accordingly.
% INPUTS:
%   rm_list (array): list of indices to remove
%   nodePos (matrix): [x,3] matrix of the (x,y,z) coordinates for
%                           each node in the graph. total number nodes = x
%   nodeEdges (matrix): [y,2] matrix. Each row corresponds to an edge
%                             in the graph, and the two values equal the
%                             node index that the edge connects. For
%                             example, if the first entry is [1,2], then
%                             this indicates the first edge is connecting
%                             node 1 and node 2.
% OUTPUTS:
%   pos_reindex (): list of node positions, excluding elments in rm_list
%   edge_reindex (): edge/node mapping. Reindexed to account for removed
%                       elments.

% Initialize 
nNodes = size(nodePos,1);
map = (1:nNodes)';
map(rm_list) = [];

% 
mapTemp = (1:length(map))';
nodeMap = zeros(nNodes,1);
nodeMap(map) = mapTemp;

% Reindex edge matrix
edgesNew = nodeMap(nodeEdges);
[ir,~] = find(edgesNew == 0);

edgesNew(ir,:) = [];
nodePos(rm_list,:) = [];

% Set output variables
edge_reindex = edgesNew;
pos_reindex = nodePos;

end

