function [pos_unique, edges_mapped] =...
downsample_loops(loop_node_idcs, nodes, edges, delta, protect)
%%regraphNodes_new Downsample a graph
% INPUTS
%   loop_node_idcs (cell array): [1,1,N] each entry contains an array of
%                            indices that can be down sampled.
%   nodes (array): [X, Y, Z] matrix of coordinates for nodes
%   edges (array): [M, 2] matrix of edges connecting node indices
%   delta (uint): step size for searching
%   protect (boolean): 
%                      True = protect end points.
%                      False = down sample end points.
%
% OUTPUTS
%   nodes_out (array): coordinates of nodes after down sampling
%   edges_out (array): matrix of edges after down sampling
%{
--------------------------------------------------------------------------
PURPOSE:
This function was designed to down sample nodes belonging to loops. It is
meant to avoid down sampling non-loop strucutres, to retain the original
vascular morphology.
--------------------------------------------------------------------------
TODO:
1) Some loops are not being downsampled, even when the delta is larger
than the distance between nodes. This may be due to some nodes remaining
protected, thus preventing complete removal of all loops. Need to
investigate futher.
--------------------------------------------------------------------------
%}

%% Re-index the edges containing the nodes to downsample
%%% Create list of node indices (from all segment groups) to downsample.
% This will be used to find all edges involved in downsampling, which will
% be used to create a separate edge list.
nodes_ds = unique(horzcat(loop_node_idcs{:}));
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

%%% Verify that all nodes in loop_node_idcs are contained in edges_ds
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

%% Regraph (downsample)
%{
CODE EXPLANATION:

Inside (for ii=2:n_nodes), this function performs the following
    For the current node in the for loop, find the redundant nodes within a
    search radius (+/- [delta, delta, delta]).
    
    if there are no nodes within radius:
        add the current node position to the array "pos_unique", which will be
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
    nodes_out = pos_unique;
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
% the non-loop nodes into "pos_unique" will ensure these nodes are considered
% "unique," and they will not be down sampled.
%
% Create array of node indices to keep
nkeep = 1:1:size(nodes,1);
nkeep = nkeep';

%%% If protecting end nodes
%       - add indices to array of nodes to keep.
%       - remove from indices of nodes to downsample.
% Otherwise nkeep will only be the nodes in loops.
if protect
    % Remove end node indices from array of node indices to down sample
    idcs = ~ismember(end_nodes, nodes_ds);
    nodes_ds = nodes_ds(idcs,:);
end

%%% Array of nodes to not down sample
nkeep = ~ismember(nkeep, nodes_ds);
% Find array indices equal to 1. These are the node indices of non-loop
% segments, which should be preserved during down sampling.
nkeep = find(nkeep == 1);
% Take unique elements of nkeep to avoid duplicates
nkeep = unique(nkeep);

%%% Initialize node_map to have an index for every node
node_map = zeros(size(nodes,1), 1);
% TODO: this logic may be incorrect. May need to map the index of nkeep to
% the respective index in node_map.
node_map(nkeep) = nkeep;
% % Set the indices of the nodes to protect. The mapped index will be equal
% % to the index in the array node_map. This ensures a one-to-one mapping for
% % the protected nodes, so they will retain their original positions.
% node_map(1:length(nkeep)) = 1:length(nkeep);

%%% Visualize the graph with protected nodes
visualize_graph(nodes, edges,'Protected Nodes (green)', nkeep);

%%% "pos_unique" tracks the index positions of unique nodes.
% The current node position in the for loop (pos_tmp) is compared to all
% nodes within pos_unique. 
%
% When a node in "nodes_ds" lacks neighbors within the search radius,
% its position is added to this list, and its index is added to node_map.
% The index "n_unique" tracks the position in the "pos_unique" matrix.
%
% The nodes belonging to the non-loop edges will be categorized as unique
% nodes, since they should remain after down sampling. This is to ensure
% the down sampled loops can be reintegrated back into the non-loop
% segments.
pos_unique = nodes(nkeep,:);
n_unique = length(nkeep);
pos_unique(n_unique,:) = nodes(n_unique,:);

%%% Matrix of loop node coordinates
pos_loop = nodes(nodes_ds,:);
% TODO:
%   - may need additional variable, similar to pos_unique. 
%   - may need to update pos_loop within the if statement "isempty(red...)

%% Issue: loop nodes are not contained within "pos_unique"
% The find function compares the position of the current node "pos_tmp" to
% the nodes in pos_unique. However, pos_unique only contains protected nodes
% (not belonging to cycles).
% TODO:
%   - compare pos_tmp to the other nodes in loops (nodes_ds)
%   - reindex the output to match the node_map

%%% Iterate over all nodes & perform down sampling
for ii=1:size(nodes_ds,1)
% for ii=2:size(nodes,1)
    % Position of node under comparison
%     pos_tmp = nodes(ii,:);
    pos_tmp = pos_loop(ii,:);
    % If the node is protected, then skip this for-loop iteration
    if ~all(ismember(pos_tmp, pos_unique))
        pause(0.1);
    end
    % Find nodes within the search radius of pos_tmp
%     redundant_nodes = find(...
%         pos_tmp(1)>=(pos_unique(:,1)-hxy) & pos_tmp(1)<=(pos_unique(:,1)+hxy) & ...
%         pos_tmp(2)>=(pos_unique(:,2)-hxy) & pos_tmp(2)<=(pos_unique(:,2)+hxy) & ...
%         pos_tmp(3)>=(pos_unique(:,3)-hz) & pos_tmp(3)<=(pos_unique(:,3)+hz) );

    redundant_nodes = find(...
        pos_tmp(1)>=(pos_loop(:,1)-hxy) & pos_tmp(1)<=(pos_loop(:,1)+hxy) & ...
        pos_tmp(2)>=(pos_loop(:,2)-hxy) & pos_tmp(2)<=(pos_loop(:,2)+hxy) & ...
        pos_tmp(3)>=(pos_loop(:,3)-hz) & pos_tmp(3)<=(pos_loop(:,3)+hz) );
    
    %%% No nodes within search radius [hxy, hxy, hz] of current node.
    % Iterate n_unique & set node index in node_map(ii) equal to n_unique.
    % This will map the new node position (after down sampling) to the
    % current unique node position.
    if isempty(redundant_nodes)
        n_unique = n_unique + 1;
        node_map(ii) = n_unique;
        pos_unique(n_unique,:) = pos_tmp;
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
                d(r) = norm(pos_tmp-pos_unique(redundant_nodes(r),:));
            end
            % Retrieve index of closest node in redundant_nodes
            [~, closest_node] = min(d);
        % Otherwise the index of closest node is the first element
        else
            closest_node = 1;
        end
        %%% Replace current node with closest loop node.
        % The array node_map is for mapping all nodes in the graph.
        % The array nodes_ds contains the indices of all loop nodes (and
        % end nodes if protect is false), in the context of the entire
        % graph. Therefore, an index of X in nodes_ds corresponds to node
        % index X in the graph.
        % The following line converts the index within redundant_nodes to
        % the index within the entire graph.
        node_map(nodes_ds(ii)) = nodes_ds(redundant_nodes(closest_node));
%         node_map(ii) = redundant_nodes(closest_node);
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

%%% Array of indices representing indices of the new node positions. 
% Create array from 1 to length of pos_unique (downsampled node positions)
node_ds_idcs = 1:size(pos_unique,1);

%%% Find unconnected nodes (not contained in edges (edges_mapped))
solo = ~ismember(node_ds_idcs, edges_mapped);
solo = find(solo == 1, 1);
fprintf('Total disjoint nodes = %d\n', length(solo))

%%% If unconnected nodes, then reindex edges to exclude them
if ~isempty(solo)
    [pos_unique, edges_mapped] = rm_disjoint_nodes(pos_unique, edges_mapped);
end

%% Visual Verification

%%% Print number of removed nodes
fprintf('Regraph reduced %d nodes to %d\n',size(nodes,1),size(pos_unique,1))

%%% Verify subgraph with figure
% Highlight end nodes in visualize_graph
% highlight_nodes = node_map(1:length(nkeep));
visualize_graph(pos_unique, edges_mapped, 'After Downsampling Graph',[]);
% xlim([160, 240]); ylim([0, 80]); zlim([10,50]); view(3);

% Other nested loops
% xlim([40, 120]); ylim([180, 260]); zlim([70,110]); view(3);

end