function [edges_rm] = rm_loop_edge(nodes, edges, sp, nsp, esp)
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
%   nsp (cell array): Sparse loop nodes. Each row contains node indices in
%                       a sparse loop.
%   esp (cell array): Sparse loop edges. Each row contains edge indices in
%                       a sparse loop. Edge indices are respective to csp.
%
% OUTPUTS
%   edges_rm (matrix): edge matrix without longest edge for each loop
%
% ISSUE: the matlab function "allcycles" is supposed to return edge indices
% belonging to cycles in a graph. However, it consistently returns the
% incorrect edge indices. There is no documentation online regarding this
% issue.
% RESOLUTION:
%   - write custom function to determine edges belonging to loops
%   - Incorporate into this function

%% Remove longest edge from one cycle until all sparse cycles are opened
cnt = 1;
while any(sp)

    %% Extract the first sparse loop nodes and edges from cell array
    % Edge indices
    cedges = esp{1,:};
    % Edge matrix [source node, target node]
    cedges = edges(cedges,:);
    
    %% Calculate distances for all nodes in graph
    % Find difference between each point
    dmat = pdist(nodes, 'euclidean');
    % Convert difference to upper triangle square matrix (less memory)
    dmat = triu(squareform(dmat),1);

    %%% Find dmat indices corresponding to edges in the cycle
    % Convert indices to subscript
    sb = sub2ind(size(dmat), cedges(:,1), cedges(:,2));

    %%% Set non-edge matrix entries equal to zero
    % Array of all indices of dmat
    del_idx = 1: (size(dmat,1) * size(dmat,2));
    % Remove indices corresponding to real edges
    del_idx(sb) = [];
    % Set non-edge distances equal to zero
    dmat(del_idx) = 0;
    
    %% Old method
    %{
    %% Calculate euclidean distance between nodes in edges
    % Extract coordinates of each node in loop
    pmat = nodes(cnodes,:);
    % Find difference between each point
    dmat = pdist(pmat, 'euclidean');
    % Convert difference to upper triangle square matrix (less memory)
    dmat = triu(squareform(dmat),1);
    
    %% Remove entries that do not correspond to edges
    % The dmat matrix contains Euclidean distances between nodes that are
    % not connected by an edge. These should be discarded before finding
    % the longest edge. The nodes in the edges are indexed according to the
    % entire graph, but the nodes in dmat are indexed to just the cycle. It
    % is necessary to reindex the nodes to dmat.

    %%% Find dmat indices corresponding to edges in the cycle
    % Ordered list of node indices in edges
    old = sort(unique(cedges(:)));
    % New list of indices (1 : N)
    new = 1:1:length(old);
    % Reindex edges
    ch = changem(cedges,new, old);
    % Convert indices to subscript
    sb = sub2ind(size(dmat), ch(:,1), ch(:,2));
    
    %%% Set non-edge matrix entries equal to zero
    % Array of all indices of dmat
    del_idx = 1: (size(dmat,1) * size(dmat,2));
    % Remove indices corresponding to real edges
    del_idx(sb) = [];
    % Set non-edge distances equal to zero
    dmat(del_idx) = 0;
    %}
   
    %% Remove long edge from "edges" matrix
    
    %%% Find end nodes of longest edge
    % Find dmat index of longest edge (longest eucl. dist. to other nodes)
    [~, edge_idx] = max(dmat,[],"all","linear");
    % Convert matrix index to subscript
    [r, c] = ind2sub(size(dmat), edge_idx);
    % Find edge containing start/end nodes
    e_idx = (edges(:,1)==r & edges(:,2)==c);
    
    %%% Remove long edge from "edges" matrix
    edges(e_idx,:) = [];

    %% Recalculate graph sparsity for while-loop
    sp = graph_sparsity(edges);
    % Generate graph
    g = graph(edges(:,1), edges(:,2));
    % Find cycles in graph
    [~, esp] = allcycles(g);
    % Convert sparsity array to boolean
    sp = boolean(sp);
    % Keep edge indices from sparse cycles
    esp(~sp) = [];
    
    %% Plot updated graph
    tstr = strcat('Iteration ', num2str(cnt));
    visualize_graph(nodes, edges, {'Longest Edge Removed',tstr},[]);
    cnt = cnt+1;
end


%%% Reassign output variable to make Matlab happy
edges_rm = edges;

end

