function [nodes, edges_reindex] =...
    remove_reindex_nodes(rm_idcs, nodes, edges)
% Remove nodes from nodes and edges. Reindex edges accordingly.
% This function is used within the GUI.
%
% INPUTS:
%   rm_idcs (array): array of node indices to remove
%   nodes (matrix): [n,3] matrix of the (x,y,z) coordinates for
%                   	each node in the graph. total number nodes = n
%   edges (matrix): [n,2] matrix. Each row corresponds to an edge
%                       in the graph, and the two values equal the
%                       node index that the edge connects. For
%                       example, if the first entry is [1,2], then
%                       this indicates the first edge is connecting
%                       node 1 and node 2.
% OUTPUTS:
%   nodes_reindex (): list of node positions, excluding elments in rm_idcs
%   edges_reindex (): edge/node mapping. Reindexed to account for removed
%                       elments.

% Initialize map
n = size(nodes,1);
map = (1:n)';
map(rm_idcs) = [];

% Initialize copy of map
mapTemp = (1:length(map))';
nodeMap = zeros(n,1);
nodeMap(map) = mapTemp;

% Reindex edge matrix
edges_reindex = nodeMap(edges);
[ir,~] = find(edges_reindex == 0);

edges_reindex(ir,:) = [];
nodes(rm_idcs,:) = [];

end

