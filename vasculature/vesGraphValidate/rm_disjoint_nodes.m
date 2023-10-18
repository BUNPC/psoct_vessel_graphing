function [nodes_re, edges_re] = rm_disjoint_nodes(nodes, edges)
%% Remove disjoint nodes (nodes without any edge connections)
% The input "nodes" contains disjoint nodes, which are not contained
% anywhere within the "edges" matrix. However, the edges matrix is indexed
% according to the "nodes" matrix, meaning that the number of indices
% exceeds the number of nodes in the graph. This raises an error when
% trying to convert the nodes/edges to the Matlab graph structure.

% This function takes a list of edges, finds all the unique node indices,
% and then reindexes the edges to exclude unconnected nodes. For example,
% consider a nodes matrix with 10 nodes, and an edges matrix that only uses
% 8 nodes:
%   nodes: [x1, y1, z1; x2, y2, z2; ...; x10, y10, z10]
%   edges: [n1, n2; n3, n4; n5, n6; n9, n10]
%
% This code will create a sorted list of node indices (n1 - n6, n9, n10).
% Then it will reindex "edges" in a manner that shifts down the upper node
% indices to create a continuous array of node indices. This excludes the
% unused nodes, leading to a new node array (nodes_re).
% 
% INPUTS:
%   nodes ((N,3) integer matrix): node positions [x,y,z]
%   edges ((M,2) integer matrix): [start node, end node]
% OUTPUTS:
%   nodes_re (integer matrix): node positions [x,y,z]
%           (N-number of disjoint nodes, 3)
%   edges_re (integer matrix): [start node, end node] reindexed to the
%           nodes_re matrix.

%%% Initialize variables for re-indexing
% Create ordered list of edge indices
edso = sort(edges(:));
% Track the new re-indexed node index
ni = 1;
% Array for storing reindexed node indices
nre = zeros(length(edso),1);
% Store previous element in for-loop for comparison to current element
prior = edso(1);
% Matrix for storing node positions (after re-indexing)
nodes_re = zeros(length(unique(edso)),3);

% Enter debug on error
dbstop if error

%%% Iterate over ordered list of node indices in edges
for e=1:length(edso)
    % If the current element of the ordered indices equals the prior elem.
    if edso(e) == prior
        % Set the reindexed node index equal to the prior node index
        nre(e) = ni;
        % Update reindexed node position
        nodes_re(ni,:) = nodes(edso(e), :);
    else
        % Increment the new reindexed node index
        ni = ni + 1;
        % Update reindexed node index
        nre(e) = ni;
        % Update reindexed node position
        nodes_re(ni,:) = nodes(edso(e), :);
    end
    % Update the value of the prior element for comparison
    prior = edso(e);
end

%%% Remap node indices in edges with the re-indexed values
% This command will replace the value edso(k) with nre(k) in edges
edges_re = changem(edges, nre, edso);

end