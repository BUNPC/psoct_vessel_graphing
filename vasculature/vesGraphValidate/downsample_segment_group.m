function [nodes_out, edges_out] =...
downsample_segment_group(group_node_idcs, nodes, edges, delta)
%%regraphNodes_new Downsample a graph
% INPUTS
%   group_idcs (cell array): [1,1,N] each entry contains an array of
%                            indices that can be down sampled.
%   nodes (array): [X, 3] matrix of coordinates for nodes
%   edges (array): [Y,2] matrix of edges connecting node indices
%   delta (uint): step size for searching
% OUTPUTS
%   nodes (array): coordinates of nodes after down sampling
%   edges (array): matrix of edges after down sampling
%
%{
TODO:
1) Add argument to downsample either:
    - individual segments to smooth single segment
    - list of segments. This is useful for regraphing a
        group of short segments, which are outputted from the
        marching ellpisoid code.
2) Re-index the nodes/edges NOT in the group_node_idcs
    - create separate arrays/matrices for storing these:
    - nodepos_og, edges_og
3) After downsampling, recombine the two groups of matrices (re-indexed &
    original).
    - 
--------------------------------------------------------------------------
CODE EXPLANATION:

The loop "for n = 1:size(group_node_idcs,3)" iterates over the group of
    node indices that meet the conditions for down sampling.

Inside (for ii=2:n_nodes), this function performs the following
    For the current node in the for loop, find the redundant nodes within a
    search radius (+/- [delta, delta, delta]).
    
    Find the redundant node that is closest to the current node in the
    loop. Reassign the index of the current node to the index of the
    closest redundant node. This is performed in the following line:

    node_map(ii) = redundant_nodes(closestNode);

    This will essentially remove the current node and replace it with the
    closest redundant node. This replacement occurs after the for loop
    completes. 

After (for ii=2:n_nodes), this function calls a line for remapping the nodes
and edges:

    edges_mapped = node_map(edges_mapped);

    The array edges_mapped is an Mx2 array containing the node indices for each
    edge. node_map is the same length as the original node array (node).
    Each entry in the array node_map corresponds to the updated index of
    each node after downsampling.

    This process results in redundant edges that connect the same nodes.
    The following section removes these redundant edges.
%}

%% Re-index the edges containing the nodes to downsample
%%% Create list of node indices (from all segment groups) to downsample.
% This will be used to find all edges involved in downsampling, which will
% be used to create a separate edge list.
nodes_ds = horzcat(group_node_idcs{:});

%%% Iterate over all nodes and find edges connected to each nodes
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

%% Re-index the nodes in group_idcs
% The node indices in group_idcs are indexed based on the non-downsampled
% graph. These indices must be reindexed to 1. Similarly, the 

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

%%% Update the node indices in edges_ds with the re-indexed values
% This command will replace the value edso(k) with nre(k) in edges_ds
edges_ds_re = changem(edges_ds, nre, edso);


%% Regraph (downsample)
% Set search delta for x,y,z
hxy = delta;
hz = delta;
% Struct for storing mapping
map = struct;

% Iterate over segment ID number. Each node has a correpsonding segment ID
for n = 1:size(group_node_idcs,3)
    % Group of nodes to downsample
    nodes_idcs = cell2mat(group_node_idcs(:,:,n));
    % Remap nodes according to reindexing in last section
    nodes_idcs = changem(nodes_idcs, nre, edso);
    % Number of nodes in group of segments
    n_nodes = length(nodes_idcs);
    
    % Initalize variables for group of segments
    % TODO: determine how to re-reference node_map
    n_unique = 1;
    node_map = zeros(n_nodes,1);
    node_map(1) = 1;
       
    % Initial value of node position for comparision
    pos_new = nodes_ds_re(nodes_idcs(1),:);

    for ii=2:n_nodes
        % Position of node under comparison
        pos_tmp = nodes(nodes_idcs(ii),:);
        % Find nodes within the search radius of pos_tmp
        redundant_nodes = find(...
            pos_tmp(1)>=(pos_new(:,1)-hxy) & pos_tmp(1)<=(pos_new(:,1)+hxy) & ...
            pos_tmp(2)>=(pos_new(:,2)-hxy) & pos_tmp(2)<=(pos_new(:,2)+hxy) & ...
            pos_tmp(3)>=(pos_new(:,3)-hz) & pos_tmp(3)<=(pos_new(:,3)+hz) );
        % No nodes within search radius [hxy, hxy, hz] of current node.
        if isempty(redundant_nodes)
            n_unique = n_unique+1;
            node_map(ii) = n_unique;
            pos_new(n_unique,:) = pos_tmp;
        % At least 1 node within search radius [hxy, hxy, hz]
        else
            % If more than one node within radius
            if length(redundant_nodes)>1
                % Initialize array to track distance between current
                % node and the nodes within search radius.
                d = zeros(length(redundant_nodes),1);
                % Iterate over list of nodes within search radius
                for r=1:length(redundant_nodes)
                    d(r) = norm(pos_tmp-pos_new(redundant_nodes(r),:));
                end
                % Retrieve index of closest node in redundant_nodes
                [~, closestNode] = min(d);
                % Delete local variable
                clear d
            % Otherwise the index of closest node is the first element
            else
                closestNode = 1;
            end
            % Replace current node with closest node
            node_map(ii) = redundant_nodes(closestNode);
        end
    end
    %%% Add mapping to struct b/c it will be overwritten on next iteration
    % Re-index node_map with the node indices
    map(n).node_map = nodes_idcs(node_map);
    map(n).pos_new = pos_new;
end

%%% Create master list of node_map for entire graph
% Convert from struct to cell array
node_map_master = struct2cell(map);
% Concatenate cell arrays into a double array
node_map_master = cell2mat(node_map_master(1,:));
% Find the new edges based upon node_map
edges_mapped = node_map_master(edges_ds_re);

%%% Create master list of new node positions
% Convert from struct to cell array
pos_new = struct2cell(map);
pos_new = pos_new(2,:);
% Concatenate cell arrays into a double array
pos_new = vertcat(pos_new{:});

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
fprintf('Regraph reduced %d nodes to %d\n',n_nodes,size(nodes_out,1))

end