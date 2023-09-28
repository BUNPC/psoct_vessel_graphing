function [nodes_out, edges_out] =...
downsample_subset(subset_node_idcs, nodes, edges, delta)
%%regraphNodes_new Downsample a graph
% INPUTS
%   subset_node_idcs (cell array): [1,1,N] each entry contains an array of
%                            indices that can be down sampled.
%   nodes (array): [X, 3] matrix of coordinates for nodes
%   edges (array): [Y,2] matrix of edges connecting node indices
%   delta (uint): step size for searching
% OUTPUTS
%   nodes (array): coordinates of nodes after down sampling
%   edges (array): matrix of edges after down sampling
%
%{
--------------------------------------------------------------------------
PURPOSE:
This function was designed to down sample nodes belonging to loops. It is
meant to avoid down sampling non-loop strucutres, to retain the original
vascular morphology.
--------------------------------------------------------------------------
TODO:
1) Re-index the nodes/edges NOT in the array subset_node_idcs
    - create separate arrays/matrices for storing these:
    - nodepos_og, edges_og
2) Identify enpoints of loops connected to non-loops
    - this is currently not finding all end point in the loops. It may
    suffice to just use the end points of the non-loop segments.
    - find edges containing one node in a loop (subset_node_idcs)
    - add these nodes to the pos_new matrix so they are categorized as
    unique.
3) After downsampling, recombine the down-sampled graph with the
    non-down-sampled graph.
    - Since we maintain the non-loop edges connected to loops, we may be
    able to reindex both sets to connect the edges.
--------------------------------------------------------------------------
%}

%% Re-index the edges containing the nodes to downsample
%%% Create list of node indices (from all segment groups) to downsample.
% This will be used to find all edges involved in downsampling, which will
% be used to create a separate edge list.
nodes_ds = unique(horzcat(subset_node_idcs{:}));
nodes_ds = nodes_ds';

%%% Iterate over all nodes and find edges connected to loops
% Variable for storing edge indices
edges_ds = [];
for n = 1:length(nodes_ds)
    % Find all edges connected to nodes_ds(n)
    idcs = find(edges(:,1)==nodes_ds(n) | edges(:,2)==nodes_ds(n));
    % Add to array for storing edge indices (edge_ds)
    edges_ds = [edges_ds; idcs];
end
% Find unique edge subscripts
edges_ds = unique(edges_ds);

%%% Create new matrix of edges (containing nodes to down sample)
% The values in edges_ds are the indices of all edges that contain at least
% one of the nodes in nodes_ds. The matrix "edges" contains all edges in
% the entire graph. This line of code will use the edge indices in edges_ds
% to create a new matrix of edges [node_start, node_end].
edges_ds = edges(edges_ds,:);

%%% Verify that all nodes in subset_node_idcs are contained in edges_ds
% Note that there will likely be nodes in edges_ds that are not contained
% within nodes_ds. A loop may contain a node, which connects via an edge to
% another node not in the loop. In this case, the non-loop node in the edge
% will be contained within "edges_ds" but not "nodes_ds"
node_diff = setdiff(unique(nodes_ds(:)), unique(edges_ds(:)));
assert(isempty(node_diff),...
    'At least one node in the node down sample list is not contained ',...
    'within the list of edges.')

%% Identify last node in segment prior to edge connected to loop.
% This preserves the edges connecting loops to non-loop segments.
% This is necessary for reintegrating the down-sampled loops back
% into the non-loop vessels.

%%% Find the end node of the non-loop segment connected to the loop.
seg_end_nodes = setdiff(unique(edges_ds(:)), unique(nodes_ds(:)));

%%% Find node(s) in loop(s) connected to the segment end-node
% There could be multiple loops connected to the same segment end point.
% "n_idcs" stores the end-node index in the loop
n_idcs = [];
for n = 1:length(seg_end_nodes)
    %%% Find array indices of all edges connected to seg_end_nodes(n)
    edge_idcs = edges_ds(:,1)==seg_end_nodes(n) | edges_ds(:,2)==seg_end_nodes(n);
    % Use edges_ds array index to retrieve node indices of end points
    node_idcs = edges_ds(edge_idcs,:);
    % Convert from 2D array to 1D vertical array to enable concatenation in
    % last step of for-loop.
    node_idcs = node_idcs(:);
    
    %%% Find end nodes belonging just to a loop
    idx = node_idcs(ismember(node_idcs, nodes_ds(:)));
    % Store node indices
    n_idcs = [n_idcs; idx];
end
% Find unique nodes
n_idcs = unique(n_idcs);
% Combine segment end nodes and loop end nodes
end_nodes = [seg_end_nodes; n_idcs];
end_nodes = unique(end_nodes);

%%% Find edges connecting non-loop segments to loop segments
% The array "end_nodes" contains the end node indices of both non-loop and
% loop segments. This finds edges where both start and end node belong to
% the set of end_nodes.
e_idcs = find(ismember(edges(:,1), end_nodes) & ismember(edges(:,2), end_nodes));
% Convert from array indices to edge values
end_edges = edges(e_idcs,:);

%%% Visualize the graph
visualize_graph(nodes, edges, 'All End Nodes (exclude from down sampling)',end_nodes);
xlim([160, 240]); ylim([0, 80]); zlim([10,50]); view(3);

%% Re-index the nodes in subset_node_idcs
% The node indices in subset_node_idcs are indexed based on the
% non-downsampled graph. These indices must be reindexed to 1. Similarly,
% the edge indices must also be reindexed to 1.
% TODO: make this into a function outside of this function. It will not be
% used here, but it may be useful later.
%{
%%% Initialize variables for re-indexing
% Create ordered list of edge indices
edso = sort(edges_ds(:));
% Track the new re-indexed node index
ni = 1;
% Array for storing reindexed node indices
nre = zeros(length(edso),1);
% Temporary variable for storing last element for comparison
last = edso(1);
% Matrix for storing node positions (after re-indexing)
nodes_ds_re = zeros(length(unique(edso)),3);

%%% Iterate over ordered list of node indices in edges_ds
for e=1:length(edso)
    % If the current element of the ordered indices equals the prior elem.
    if edso(e) == last
        % Set the reindexed node index equal to the last node index
        nre(e) = ni;
        % Update reindexed node position
        nodes_ds_re(ni,:) = nodes(edso(e), :);
    else
        % Increment the new reindexed node index
        ni = ni + 1;
        % Update reindexed node index
        nre(e) = ni;
        % Update reindexed node position
        nodes_ds_re(ni,:) = nodes(edso(e), :);
    end
    % Update the value of the last element for comparison
    last = edso(e);
end

%%% Remap node indices in edges_ds with the re-indexed values
% This command will replace the value edso(k) with nre(k) in edges_ds
edges_ds_re = changem(edges_ds, nre, edso);

%%% Remap node indices belonging to segment and loop end points
end_nodes_re = changem(end_nodes, nre, edso);

%%% Verify subgraph visually (figure) and programatically
visualize_graph(nodes_ds_re, edges_ds_re, 'Subset of graph',end_nodes_re);
xlim([160, 240]); ylim([0, 80]); zlim([10,50]); view(3);
%}
%% Regraph (downsample)
%{
CODE EXPLANATION:

Inside (for ii=2:n_nodes), this function performs the following
    For the current node in the for loop, find the redundant nodes within a
    search radius (+/- [delta, delta, delta]).
    
    if there are no nodes within radius:
        add the current node position to the array "pos_new", which will be
        the new node positions after completing the for loop.
    else:

        Find the redundant node that is closest to the current node in the
        loop. Reassign the index of the current node to the index of the
        closest redundant node. This is performed in the following line:
    
            node_map(ii) = redundant_nodes(closestNode);
    
        This will essentially remove the current node and replace it with
        the closest redundant node. This replacement occurs after the for
        loop completes. 

After (for ii=2:n_nodes), this function calls a line for remapping the nodes
and edges:

    edges_mapped = node_map(edges_mapped);
    nodes_out = pos_new;
    edges_out = edges_mapped;

    The array edges_mapped is an Mx2 array containing the node indices for each
    edge. node_map is the same length as the original node array (node).
    Each entry in the array node_map corresponds to the updated index of
    each node after downsampling.

    This process results in redundant edges that connect the same nodes.
    The following section removes these redundant edges.
%}

%%% Set search delta for x,y,z
hxy = delta;
hz = delta;

%%% Create array of nodes that will NOT be down sampled
% The goal is to preserve the morphology of the non-loop segments. Placing
% the non-loop nodes into "pos_new" will ensure these nodes are considered
% "unique," and they will not be down sampled.
%
% Create array of node indices 
nkeep = 1:1:size(nodes,1);
nkeep = nkeep';
% Find nodes not in loops
nkeep = ~ismember(nkeep, nodes_ds);
% Find array indices equal to 1. These are the node indices of non-loop
% segments, which should be preserved during down sampling.
nkeep = find(nkeep == 1);

% Add the end node indices to the array of nodes to keep
nkeep = [nkeep; end_nodes];
% Take unique elements of nkeep to avoid duplicates
nkeep = unique(nkeep);

node_map = zeros(size(nodes,1), 1);
node_map(1:length(nkeep)) = 1:length(nkeep);

%%% "pos_new" tracks the index positions of unique nodes.
% The current node position in the for loop (pos_tmp) is compared to all
% nodes within pos_new. 
%
% When a node in "nodes_ds_re" lacks neighbors within the search radius,
% its position is added to this list, and its index is added to node_map.
% The index "n_unique" tracks the position in the "pos_new" matrix.
%
% The nodes belonging to the non-loop edges will be categorized as unique
% nodes, since they should remain after down sampling. This is to ensure
% the down sampled loops can be reintegrated back into the non-loop
% segments.
pos_new = nodes(nkeep,:);
n_unique = size(pos_new, 1);
pos_new(n_unique,:) = nodes(n_unique,:);

%%% Code from reindexing method
%{
% Number of nodes in the down sampling subset
n_nodes = size(nodes_ds_re, 1);

% Struct for storing mapping
node_map = zeros(n_nodes,1);
node_map(1) = 1;

% Assign first node as unique
pos_new = nodes_ds_re(end_nodes_re,:);
n_unique = size(pos_new, 1);
pos_new(n_unique,:) = nodes_ds_re(1,:);
%}

%%% Iterate over all nodes & perform down sampling
for ii=2:size(nodes,1)
    % Position of node under comparison
    pos_tmp = nodes(ii,:);
    % Find nodes within the search radius of pos_tmp
    redundant_nodes = find(...
        pos_tmp(1)>=(pos_new(:,1)-hxy) & pos_tmp(1)<=(pos_new(:,1)+hxy) & ...
        pos_tmp(2)>=(pos_new(:,2)-hxy) & pos_tmp(2)<=(pos_new(:,2)+hxy) & ...
        pos_tmp(3)>=(pos_new(:,3)-hz) & pos_tmp(3)<=(pos_new(:,3)+hz) );
    
    %%% No nodes within search radius [hxy, hxy, hz] of current node.
    % Iterate n_unique & set node index in node_map(ii) equal to n_unique.
    % This will map the new node position (after down sampling) to the
    % current unique node position.
    if isempty(redundant_nodes)
        n_unique = n_unique + 1;
        node_map(ii) = n_unique;
        pos_new(n_unique,:) = pos_tmp;
    % At least 1 node within search radius [hxy, hxy, hz]
    else
        % If more than one node within radius
        if length(redundant_nodes)>1
            % Delete local variable
            clear d
            % Initialize array to track distance between current
            % node and the nodes within search radius.
            d = zeros(length(redundant_nodes),1);
            % Iterate over list of nodes within search radius
            for r=1:length(redundant_nodes)
                d(r) = norm(pos_tmp-pos_new(redundant_nodes(r),:));
            end
            % Retrieve index of closest node in redundant_nodes
            [~, closest_node] = min(d);
        % Otherwise the index of closest node is the first element
        else
            closest_node = 1;
        end
        % Replace current node with closest node
        node_map(ii) = redundant_nodes(closest_node);
    end
end

% Find the new edges based upon node_map
edges_mapped = node_map(edges);

%% Remove single-node edges and redundant edges
% Remove single-node edges that are connecting the same node
edges_mapped = edges_mapped(edges_mapped(:,1)~=edges_mapped(:,2),:);

% Find unique edges
sE = cell(size(edges_mapped,1),1);
for ii=1:length(edges_mapped)
    if edges_mapped(ii,1)<edges_mapped(ii,2)
        sE{ii} = sprintf('%05d%05d',edges_mapped(ii,1),edges_mapped(ii,2));
    else
        sE{ii} = sprintf('%05d%05d',edges_mapped(ii,2),edges_mapped(ii,1));
    end
end

% Find indices of unique edges
[~,unique_edge_idx,~] = unique(sE);
% Remove redundant edges
edges_mapped = edges_mapped(sort(unique_edge_idx),:);

%% check for new dangling nodes (i.e. nodes with nB=1 that were nB=2
% if we want to implement this, we just need to have nB_old and nB_new and
% use the node_map to find when nB_new=1 and nB_old=2 for given nodes and
% then delete all nodes and edges back to the bifurcation node

%% Reassign output variables
nodes_out = pos_new;
edges_out = edges_mapped;
fprintf('Regraph reduced %d nodes to %d\n',size(nodes,1),size(nodes_out,1))
%%% Verify subgraph with figure
visualize_graph(nodes_out, edges_out, 'After Downsampling Graph',[]);
xlim([160, 240]); ylim([0, 80]); zlim([10,50]); view(3);
pause(0.01)

end