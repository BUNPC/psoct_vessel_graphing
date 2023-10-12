function [edges_rm] = rm_loop_edge(loop_node_idcs, nodes, edges)
%rm_loop_edge Remove longest edge of loop.
% -------------------------------------------------------------------------
% PURPOSE:
% This function handles the case where there are more than two segments
% involved in a loop. In this case, downsampling will just replace each
% endpoint with another node further away from the loop, thus expanding the
% loop. In this case, the longest edge will be removed, removing the loop.
% -------------------------------------------------------------------------
% TODO:
% - Determine if all end point nodes should be connected to one end node.
% In this scenario, one of the nodes must be designated as the "anchor" and
% the other end points must be connected to it. Then, the edges without
% connections to the end points must be removed. This would remove some of
% the morphology of the branching, so it must be compared to the
% segmentation.
% -------------------------------------------------------------------------
% INPUTS
%   loop_node_idcs (cell array): [1,1,N] each entry contains an array of
%                            indices that can be down sampled.
%   nodes (array): [X, Y, Z] matrix of coordinates for nodes
%   edges (array): [M, 2] matrix of edges connecting node indices
%
% OUTPUTS
%   edges_rm (matrix): edge matrix without longest edge for each loop
%

for ii = 1:size(loop_node_idcs, 1)
    % Find end points for each loop
    nidx_loop = loop_node_idcs{ii};
    
    %%% Calculate euclidean distance between nodes in loop
    % Extract coordinates of each node in loop
    pmat = nodes(nidx_loop,:);
    % Find difference between each point
    dmat = pdist(pmat, 'euclidean');
    % Convert difference to lower triangle square matrix
    dmat = tril(squareform(dmat),-1);

    %%% Find end nodes of longest edge
    % Find longest edge (longest eucl. dist. to other nodes)
    [~, edge_idx] = max(dmat,[],"all","linear");
    % Convert edge index to dmat subscripts
    [r, c] = ind2sub(size(dmat), edge_idx);
    % Convert dmat subscripts to node indices in "nodes"
    e1_idx = nidx_loop(r);
    e2_idx = nidx_loop(c);

    %%% Remove long edge from "edges" matrix
    % Find long edge in "edges" matrix
    e_idx = ((edges(:,1)==e1_idx & edges(:,2)==e2_idx) |...
        (edges(:,1)==e2_idx & edges(:,2)==e1_idx));
    % Remove edge from "edges" matrix
    edges(e_idx,:) = [];
       
    %%% TODO:
    %   - Test removing other edges in loop and connecting end points to anchor.
    % Steps:
    %   - Define anchor point
    %   - Remove loop edges NOT connected to anchor
    %   - Reindex edges and nodes

end

%%% Verify that loops are closed
% Create graph
g = graph(edges(:,1), edges(:,2));
% Find cycles
[~, edgecycles] = allcycles(g);
% Assert that no cycles remain
assert(isempty(edgecycles), 'Loops remain after removing longest edge.')

%%% Reassign output variable to make Matlab happy
edges_rm = edges;

end