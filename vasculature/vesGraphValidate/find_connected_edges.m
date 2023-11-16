function [edges_conn] = find_connected_edges(loop_nodes, edges)
%find_connected_edges find all edges in and connected to loop.
% Iterate over node indices (belonging to a loop)
%   - find all edges connceted to these node indices
%   - create a subgraph from edges and nodes
%       - identify which edges belong to loop
%       - identify which edges are attached to loop
% INPUTS:
%       loop_nodes (array): node indices belonging to loop
%       edges (Nx2 array): array of edges [n1, n2]
% OUTPUTS:
%       edges_conn (array): indices of all edges connected to the nodes in
%           the loop_nodes array.

%% Iterate over all nodes and find edges connected to loops
% Variable for storing edge indices
edges_conn = [];

% Iterate over node indices in loop
for n = 1:length(loop_nodes)
    % Find all edges connected to nidx_ds(n)
    idcs = find(edges(:,1)==loop_nodes(n) | edges(:,2)==loop_nodes(n));
    % Add to array for storing edge indices (edge_ds)
    edges_conn = [edges_conn; idcs]; %#ok<AGROW> 
end

% The above logic may result in repeats of the same index. This line will
% remove repeated indices.
edges_conn = unique(edges_conn);

end