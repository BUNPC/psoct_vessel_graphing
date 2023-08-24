function [group_idcs, nodes, edges, validated_nodes,validatedEdges] =...
regraphNodes_new(group_idcs, nodes, edges, validated_nodes, delta)
%%regraphNodes_new Downsample a graph
% INPUTS
%   group_idcs (cell array): [1,1,N] each entry contains an array of
%                            indices that can be merged.
%   nodes (array): [X, 3] matrix of coordinates for nodes
%   edges (array): [Y,2] matrix of edges connecting node indices
%   validated_nodes (array): binary array (1=validated)
%   delta (uint): step size for searching
% OUTPUTS
%   nodes (array): coordinates of nodes after down sampling
%   edges (array): matrix of edges after down sampling
%   validated_nodes (array):
%   validatedEdges (array):
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

    nodeEdges = node_map(nodeEdges);

    The array nodeEdges is an Mx2 array containing the node indices for each
    edge. node_map is the same length as the original node array (node).
    Each entry in the array node_map corresponds to the updated index of
    each node after downsampling.

    This process results in redundant edges that connect the same nodes.
    The following section removes these redundant edges.
%}

%% Initialization
% Reassign nodes/edges to local variables
nodeEdges = edges;

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
    
    % Subset of validated nodes for group node indices
    group_validated_nodes = validated_nodes(nodes_idcs);
    validated_nodes_new(1) = group_validated_nodes(1);
       
    % Initial value of node position for comparision
    pos_new = nodes(nodes_idcs(1),:);

    for ii=2:n_nodes
        % Position of node under comparison
        pos_tmp = nodes(nodes_idcs(ii),:);

        % If this node has not been validated
        if validated_nodes(ii)==0  
            redundant_nodes = find(...
                pos_tmp(1)>=(pos_new(:,1)-hxy) & pos_tmp(1)<=(pos_new(:,1)+hxy) & ...
                pos_tmp(2)>=(pos_new(:,2)-hxy) & pos_tmp(2)<=(pos_new(:,2)+hxy) & ...
                pos_tmp(3)>=(pos_new(:,3)-hz) & pos_tmp(3)<=(pos_new(:,3)+hz) );
            
            %%% This is old code which doesn't work.
            %{
            % in the redundant_nodes remove unconnected nodes to current processing node. This will
            % avoid unwanted connections and loops
            temp_lst = redundant_nodes;
            new_lst = [];
            curr_node = ii;
            for u = 1:length(redundant_nodes)
               node_edges = edges(:,1) == curr_node | edges(:,2) == curr_node;
               conn_nodes = edges(node_edges,:);
               conn_nodes = conn_nodes(:);
               common_node = intersect(temp_lst,conn_nodes);
               if isempty(common_node)
                   break;
               else
                  new_lst = [new_lst; common_node];
                  curr_node = common_node(1);
                  temp_lst = setdiff(temp_lst,curr_node);
               end
            end
            redundant_nodes = new_lst;
            %}
            %%%
            
            % No nodes within search radius [hxy, hxy, hz] of current node.
            % Updated the new list of validated nodes with a 0 for this
            % node ID.
            if isempty(redundant_nodes)
                valid_flag = 0;
                n_unique = n_unique+1;
                node_map(ii) = n_unique;
                pos_new(n_unique,:) = pos_tmp;
                validated_nodes_new(n_unique) = valid_flag;
            
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
       
        % If this is a validated node
        else
            valid_flag = 1;
            n_unique = n_unique + 1;
            node_map(ii) = n_unique;
            pos_new(n_unique,:) = pos_tmp;
            validated_nodes_new(n_unique) = valid_flag;
        end
        %%% Create a copy of the node_map, 
        map(n).node_map = node_map;
        map(n).pos_new = pos_new;
        map(n).validated_nodes_new = validated_nodes_new;
    end
end

%%% Find the new edges based upon node_map
nodeEdges = node_map(nodeEdges);

%% Remove single-node edges and redundant edges
% Remove single-node edges that are connecting the same node
nodeEdges = nodeEdges(nodeEdges(:,1)~=nodeEdges(:,2),:);

% Find unique edges
sE = cell(size(nodeEdges,1),1);
for ii=1:length(nodeEdges)
    if nodeEdges(ii,1)<nodeEdges(ii,2)
        sE{ii} = sprintf('%05d%05d',nodeEdges(ii,1),nodeEdges(ii,2));
    else
        sE{ii} = sprintf('%05d%05d',nodeEdges(ii,2),nodeEdges(ii,1));
    end
end

% Find indices of unique edges
[~,unique_edge_idx,~] = unique(sE);
% Remove redundant edges
nodeEdges = nodeEdges(sort(unique_edge_idx),:);

%% check for new dangling nodes (i.e. nodes with nB=1 that were nB=2
% if we want to implement this, we just need to have nB_old and nB_new and
% use the node_map to find when nB_new=1 and nB_old=2 for given nodes and
% then delete all nodes and edges back to the bifurcation node

%% Reassign output variables
nodes = pos_new;
validated_nodes = validated_nodes_new';
edges = nodeEdges;
fprintf('Regraph reduced %d nodes to %d\n',n_nodes,size(nodes,1))

%% Assign the validated edges index
validatedEdges = zeros(size(edges,1),1);
for uu = 1:size(edges,1)
    if validated_nodes(edges(uu,1)) == 1 && validated_nodes(edges(uu,2)) == 1
        validatedEdges(uu) = 1;
    end
end

%% Update node map for validated/unvalidated nodes
%{
    function [n_unique, node_map, pos_new, validated_nodes_new, node_unique] =...
        remap(n_unique, node_map, pos_new, validated_nodes_new, ii, pos, node_unique)
n_unique = n_unique+1;
node_map(ii) = n_unique;
pos_new(n_unique,:) = pos;
validated_nodes_new(n_unique) = 0;
end
%}
end