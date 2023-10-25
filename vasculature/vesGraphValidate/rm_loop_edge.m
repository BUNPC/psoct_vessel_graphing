function [edges_rm] = rm_loop_edge(nodes, edges, sp, csp)
%rm_loop_edge Remove longest edge of loop.
% -------------------------------------------------------------------------
% PURPOSE:
% This function handles the case where there are more than two segments
% involved in a loop. In this case, downsampling will just replace each
% endpoint with another node further away from the loop, thus expanding the
% loop. In this case, the longest edge will be removed, removing the loop.
% -------------------------------------------------------------------------
% INPUTS
%   g (graph): the original graph struct
%   nodes (array): [X, Y, Z] matrix of coordinates for nodes
%   edges (array): [M, 2] matrix of edges connecting node indices
%   csp (cell array): Each row in the cell array contains the node
%           indices for a sparse loop.
%
% OUTPUTS
%   edges_rm (matrix): edge matrix without longest edge for each loop

%% TODO:
% In the case of adjoined sparse cycles, a for-loop may result in removing
% more edges than intended. Instead, a while loop will ensure that the
% sparsity matrix is regenerated between iterations and each loop is
% handled individually. 
%
%   - Change to while loop:
%       - remove longest edge of loop
%       - call graph_sparsity again
%   - Replace logic with builtin "distances" function (faster).
%   - Test removing other edges in loop and connecting end points to anchor.
%       - Define anchor point
%       - Remove loop edges NOT connected to anchor
%       - Reindex edges and nodes

%% Remove longest edge in each loop

% While at least one cycle is still sparse
while any(sp)
    % Extract the first sparse loop from cell array
    cnodes = csp{1,:};
   
    %%% Calculate euclidean distance between nodes in edges
    % TODO: this is calculating distance between all nodes (even
    % unconnected ones). Need to modify to only calculate distance between
    % nodes defined in edges.

    % Extract coordinates of each node in loop
    pmat = nodes(cnodes,:);
    % Find difference between each point
    dmat = pdist(pmat, 'euclidean');
    % Convert difference to lower triangle square matrix
    dmat = tril(squareform(dmat),-1);

    %%% Find end nodes of longest edge
    % Find dmat index of longest edge (longest eucl. dist. to other nodes)
    [~, edge_idx] = max(dmat,[],"all","linear");
    % Convert matrix index to subscript
    [r, c] = ind2sub(size(dmat), edge_idx);
    % Convert from matrix subscript to node index
    e1 = cnodes(r);
    e2 = cnodes(c);

    %%% Remove long edge from "edges" matrix
    % Find long edge in "edges" matrix
    e_idx = ((edges(:,1)==e1 & edges(:,2)==e2) |...
        (edges(:,1)==e2 & edges(:,2)==e1));
    % Remove edge from "edges" matrix
    edges(e_idx,:) = [];

end


%%% Reassign output variable to make Matlab happy
edges_rm = edges;

end