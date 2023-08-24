function [nodes, edges] =...
downsample_segment_group(group_idcs, nodes, edges, delta)
%%regraphNodes_new Downsample a graph
% INPUTS
%   group_idcs (cell array): [1,1,N] each entry contains an array of
%                            indices that can be merged.
%   nodes (array): [X, 3] matrix of coordinates for nodes
%   edges (array): [Y,2] matrix of edges connecting node indices
%   delta (uint): step size for searching
% OUTPUTS
%   nodes (array): coordinates of nodes after down sampling
%   edges (array): matrix of edges after down sampling
%
% TODO:
%   1) Iterate over list of nodes. This will be used to regraph:
%           - individual segments. Useful for smoothing single segment.
%           - list of segments. This is useful for regraphing a
%               group of short segments, which are outputted from the
%               marching ellpisoid code.
%{
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

%% Initialization
% Set search delta for x,y,z
hxy = delta;
hz = delta;

% Struct for storing mapping
map = struct;

%% Regraph (downsample)

% Iterate over segment ID number. Each node has a correpsonding segment ID
for n = 1:size(group_idcs,3)
    % Group of nodes to downsample
    nodes_idcs = cell2mat(group_idcs(:,:,n));
    % Number of nodes in group of segments
    n_nodes = length(nodes_idcs);
    
    % Initalize variables for group of segments
    % TODO: determine how to re-reference node_map
    n_unique = 1;
    node_map = zeros(n_nodes,1);
    node_map(1) = 1;
       
    % Initial value of node position for comparision
    pos_new = nodes(nodes_idcs(1),:);

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
edges_mapped = node_map_master(edges);

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
nodes = pos_new;
edges = edges_mapped;
fprintf('Regraph reduced %d nodes to %d\n',n_nodes,size(nodes,1))

end