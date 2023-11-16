function [nodes_keep_re, edges_mapped_re] =...
downsample_loops(loop_node_idcs, nodes, edges, delta, protect)
%%regraphNodes_new Downsample a graph
% INPUTS
%   loop_node_idcs (cell array): [1,1,N] each entry contains an array of
%                            indices that can be down sampled.
%   nodes (array): [X, Y, Z] matrix of coordinates for nodes
%   edges (array): [M, 2] matrix of edges connecting node indices
%   delta (uint): step size for searching
%   protect (boolean): 
%           True = protect end points.
%           False = down sample end points.
%   loop_close (boolean):
%           True = connect end points of loops
%           False = normal down sampling of loops
%
% OUTPUTS
%   nodes_keep_re (array): coordinates of nodes after down sampling and
%           reindexing to exclude hanging nodes
%   edges_mapped_re (array): matrix of edges after down sampling and
%           reindexing to remove hanging nodes
%{
--------------------------------------------------------------------------
PURPOSE:
This function was designed to down sample nodes belonging to loops. It is
meant to avoid down sampling non-loop strucutres, to retain the original
vascular morphology.
--------------------------------------------------------------------------
%}

%% Re-index the edges containing the nodes to downsample
%%% Create list of node indices (from all segment groups) to downsample.
% This will be used to find all edges involved in downsampling, which will
% be used to create a separate edge list.

% If loop_node_idcs is the indices for a single loop
if isa(loop_node_idcs, 'double')
    nidx_ds = loop_node_idcs';
% Otherwise, it is a cell array for all loops
else
    nidx_ds = unique(horzcat(loop_node_idcs{:}));
    nidx_ds = nidx_ds';
end

%%% Iterate over all nodes and find edges connected to loops
edges_ds = find_connected_edges(nidx_ds, edges);

%%% Create new matrix of edges (containing nodes to down sample)
% The values in edges_ds are the indices of all edges that contain at least
% one of the nodes in nidx_ds. The matrix "edges" contains all edges in
% the entire graph. This line of code will use the edge indices in edges_ds
% to create a new matrix of edges [node_start, node_end].
edges_ds = edges(edges_ds,:);

%%% Verify that all nodes in loop_node_idcs are contained in edges_ds
% Note that there will likely be nodes in edges_ds that are not contained
% within nidx_ds. A loop may contain a node, which connects via an edge to
% another node not in the loop. In this case, the non-loop node in the edge
% will be contained within "edges_ds" but not "nidx_ds"
node_diff = setdiff(unique(nidx_ds(:)), unique(edges_ds(:)));
assert(isempty(node_diff),...
    'At least one node in the node down sample list is not contained ',...
    'within the list of edges.')

%% Identify last node in segment prior to edge connected to loop.
% This preserves the edges connecting loops to non-loop segments.
% This is necessary for reintegrating the down-sampled loops back
% into the non-loop vessels.

%%% Find the end node of the non-loop segment connected to the loop.
seg_end_nodes = setdiff(unique(edges_ds(:)), unique(nidx_ds(:)));

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
    idx = node_idcs(ismember(node_idcs, nidx_ds(:)));
    % Store node indices
    n_idcs = [n_idcs; idx];
end
% Find unique nodes
n_idcs = unique(n_idcs);
% Combine segment end nodes and loop end nodes
end_nodes = [seg_end_nodes; n_idcs];
end_nodes = unique(end_nodes);

%%% Visualize the graph with protected nodes
% visualize_graph(nodes, edges,'Protected Nodes (green)', end_nodes);

%% Setup arrays for down sampling
%{
CODE EXPLANATION:

Inside (for ii=2:n_nodes), this function performs the following
    For the current node in the for loop, find the redundant nodes within a
    search radius (+/- [delta, delta, delta]).
    
    if there are no nodes within radius:
        add the current node position to the array "nodes_keep", which will be
        the new node positions after completing the for loop.
    else:

        Find the redundant node that is closest to the current node in the
        loop. Reassign the index of the current node to the index of the
        closest redundant node. This is performed in the following line:
    
            node_map(ii) = nkeep_close(closestNode);
    
        This will essentially remove the current node and replace it with
        the closest redundant node. This replacement occurs after the for
        loop completes. 

After (for ii=2:n_nodes), this function calls a line for remapping the nodes
and edges:

    edges_mapped = node_map(edges_mapped);
    nodes_out = nodes_keep;
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
% the non-loop nodes into "nodes_keep" will ensure these nodes are considered
% "unique," and they will not be down sampled.

% Create array of node indices to keep
nidx_keep = 1:1:size(nodes,1);
nidx_keep = nidx_keep';

%%% If protecting end nodes
%       - add indices to array of nodes to keep.
%       - remove from indices of nodes to downsample.
% Otherwise nidx_keep will only be the nodes in loops.
if protect
    % Remove end node indices from array of node indices to down sample
    nidx_ds = setdiff(nidx_ds, end_nodes);
    % Verify that end nodes were removed
    assert(~all(ismember(end_nodes, nidx_ds)),...
        'Not all end nodes removed from down sample list.');
end

%%% Array of node indices to preserve when down sampling
% Remove down sample node indices from list of nodes to keep
nidx_keep = setdiff(nidx_keep, nidx_ds);

%%% Initialize node_map to have an index for every node
node_map = zeros(size(nodes,1), 1);
% Set nidx_keep indices equal to 1
node_map(nidx_keep) = nidx_keep;

%%% Visualize the graph with protected nodes
% visualize_graph(nodes, edges,'Protected Nodes (green)', nidx_keep);
% xlim([270, 300]); ylim([140, 170]); zlim([0,5]); view(3);

%%% "nodes_keep" tracks the index positions of unique nodes.
% The current node position in the for loop (pos_tmp) is compared to all
% nodes within nodes_keep. 
%
% The nodes belonging to the non-loop edges will be categorized as unique
% nodes, since they should remain after down sampling. This is to ensure
% the down sampled loops can be reintegrated back into the non-loop
% segments.
nodes_keep = nodes(nidx_keep,:);

%% Down Sample Loops
%%% Iterate over all nodes & perform down sampling
for ii=1:size(nidx_ds,1)
% for ii=2:size(nodes,1)
    % Position of node under comparison
    pos_tmp = nodes(nidx_ds(ii),:);
    % If the node is protected, then skip this for-loop iteration
    if ~all(ismember(pos_tmp, nodes_keep))
        pause(0.1);
    end
    % Find nodes within the search radius of pos_tmp
    nkeep_close = find(...
        pos_tmp(1)>=(nodes_keep(:,1)-hxy) & pos_tmp(1)<=(nodes_keep(:,1)+hxy) & ...
        pos_tmp(2)>=(nodes_keep(:,2)-hxy) & pos_tmp(2)<=(nodes_keep(:,2)+hxy) & ...
        pos_tmp(3)>=(nodes_keep(:,3)-hz) & pos_tmp(3)<=(nodes_keep(:,3)+hz) );    
    %%% No nodes within search radius [hxy, hxy, hz] of current node.
    % Categorize current node as "unique" and place into nodes_keep
    if isempty(nkeep_close)
        % Set node_map position to its current index (one-to-one mapping)
        node_map(nidx_ds(ii)) = nidx_ds(ii);
        % Add loop node position to list of nodes to keep
        nodes_keep(end + 1,:) = pos_tmp;
    % At least 1 node within search radius [hxy, hxy, hz]
    else
        % If more than one node within radius
        if length(nkeep_close)>1
            % Delete local variable
            clear d
            % Initialize array to track distance between current
            % node and the nodes within search radius.
            d = zeros(length(nkeep_close),1);
            % Iterate over list of nodes within search radius
            for r=1:length(nkeep_close)
                d(r) = norm(pos_tmp-nodes_keep(nkeep_close(r),:));
            end
            % Retrieve index of closest node in nkeep_close
            [~, closest_node] = min(d);
        % Otherwise the index of closest node is the first element
        else
            closest_node = 1;
        end
        %%% Replace current node with closest loop node (prior version)
        % The array node_map is for mapping all nodes in the graph.
        % The array nidx_ds contains the indices of all loop nodes (and
        % end nodes if protect is false), in the context of the entire
        % graph. Therefore, an index of X in nidx_ds corresponds to node
        % index X in the graph.
        % The following line converts the index within nkeep_close to
        % the index within the entire graph.
%         node_map(nidx_ds(ii)) = nidx_ds(nkeep_close(closest_node));
        
        %%% Replace current node with closest loop node (newer version)
        % The index in nkeep_close(closest_node) corresponds to the
        % index in "nidx_keep", which is set by the line: nodes_keep = nodes(nidx_keep,:)
        % Therefore, the index "nidx_keep" must be remapped back to the
        % original index in "nodes".
        %
        % Position of closest node
        pos = nodes_keep(nkeep_close(closest_node), :);
        % Find position within entire list of nodes
        [~, nidx] = ismember(pos, nodes, 'rows');
        % Add index original node index to node map
        node_map(nidx_ds(ii)) = nidx;

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

%% Check for unconnected nodes

%%% Array of indices representing indices of the mapped positions. 
idcs_re = 1:max(node_map(:));

%%% Find hanging nodes (not contained in edges_mapped)
solo = setdiff(idcs_re, edges_mapped(:));
fprintf('Total disjoint nodes = %d\n', length(solo))

%%% If unconnected nodes, then reindex edges to exclude them
if ~isempty(solo)
    [nodes_keep_re, edges_mapped_re] = rm_disjoint_nodes(nodes, edges_mapped);
    %%% Verify the edges were reindexed 
    % Array of indices for remapped nodes
    nre = 1:size(nodes_keep_re,1);
    assert(isempty(setdiff(nre, edges_mapped_re)),...
        'Nodes and edges incorrectly reindexed.');
else
    %%% Assign output variable
    nodes_keep_re = nodes_keep;
    edges_mapped_re = edges_mapped;
end

%%% Print number of removed nodes
fprintf('Regraph reduced %d nodes to %d\n',size(nodes,1),size(nodes_keep_re,1))

pause(0.001)
end