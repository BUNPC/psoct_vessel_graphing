function [nodes, edges, validatedNodes,validatedEdges] =...
regraphNodes_new(nodes, edges,validatedNodes, delta)
%%regraphNodes_new Downsample a graph
% INPUTS
%   segn (array): The segment index of the respective node.
%   nodes (array):
%   edges (array):
%   validatedNodes (array):
%   delta (uint):
% OUTPUTS
%   nodes (array):
%   edges (array):
%   validatedNodes (array):
%   validatedEdges (array):
%
% TODO:
% - iterate over edges and only downsample nodes within edge
%       - may need to reindex nodeEdges
% - recalculate diameter (this script does not update diameters)
%{

Inside (for ii=2:nNodes), this function performs the following
    For the current node in the for loop, find the redundant nodes within a
    search radius (+/- [delta, delta, delta]).
    
    Find the redundant node that is closest to the current node in the
    loop. Reassign the index of the current node to the index of the
    closest redundant node. This is performed in the following line:

    nodeMap(ii) = redundant_nodes(closestNode);

    This will essentially remove the current node and replace it with the
    closest redundant node. This replacement occurs after the for loop
    completes. 

After (for ii=2:nNodes), this function calls a line for remapping the nodes
and edges:

    nodeEdges = nodeMap(nodeEdges);

    The array nodeEdges is an Mx2 array containing the node indices for each
    edge. nodeMap is the same length as the original node array (node).
    Each entry in the array nodeMap corresponds to the updated index of
    each node after downsampling.

    This process results in redundant edges that connect the same nodes.
    The following section removes these redundant edges.
%}

%% Initialization
% Reassign nodes/edges to local variables
nodePos = nodes;
nodeEdges = edges;
% Count number of nodes and edges
nNodes = size(nodePos,1);
nEdges = size(nodeEdges,1);
if ~exist('nodeDiam')
    nodeDiam = zeros(nNodes,1);
end

% Set search delta for x,y,z
hxy = delta;
hz = delta;

% Initalize counter of unique nodes, node map, and list of unique nodes
nNodesUnique = 1;
nodeMap = zeros(nNodes,1);
nodeUnique = zeros(nNodes,1);

% Set inital values of node map, node position, and unique nodes
nodeMap(1) = 1;
nodePosNew = nodePos(1,:);
nodeUnique(1) = 1;
validatedNodesNew(1) = validatedNodes(1);

%% Regraph (downsample)

hwait = waitbar(0,'Downsampling nodes...');   
for ii=2:nNodes
    % Waitbar update
    if isequal(rem(ii,1000),0)
        waitbar(ii/nNodes,hwait);   %updateing waitbar takes a long time. 
    end
   
    % Temporary variable for node position
    pos = nodePos(ii,:);

    % If this node has not been validated
    if validatedNodes(ii)==0  
        redundant_nodes = find(pos(1)>=(nodePosNew(:,1)-hxy) & pos(1)<=(nodePosNew(:,1)+hxy) & ...
            pos(2)>=(nodePosNew(:,2)-hxy) & pos(2)<=(nodePosNew(:,2)+hxy) & ...
            pos(3)>=(nodePosNew(:,3)-hz) & pos(3)<=(nodePosNew(:,3)+hz) );
        
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
        
        % No nodes within search radius [hxy, hxy, hz] of current node
        if isempty(redundant_nodes)
            nNodesUnique = nNodesUnique+1;
            nodeMap(ii) = nNodesUnique;
            nodeUnique(ii) = 1;
            nodePosNew(nNodesUnique,:) = pos;
            nodeDiamNew(nNodesUnique) = nodeDiam(ii);
            validatedNodesNew(nNodesUnique) = 0;
        
        % At least 1 node within search radius [hxy, hxy, hz]
        else
            % If more than one node within radius
            if length(redundant_nodes)>1
                % Initialize array to track distance between current
                % node and the nodes within search radius.
                d = zeros(length(redundant_nodes),1);
                % Iterate over list of nodes within search radius
                for r=1:length(redundant_nodes)
                    d(r) = norm(pos-nodePosNew(redundant_nodes(r),:));
                end
                % Retrieve index of closest node in redundant_nodes
                [~, closestNode] = min(d);
                % Delete local variable
                clear d
            else
                closestNode = 1;
            end
            % Replace current node with closest node
            nodeMap(ii) = redundant_nodes(closestNode);
        end
   
    % If this is a validated node
    else
        nNodesUnique = nNodesUnique+1;
        nodeMap(ii) = nNodesUnique;
        nodeUnique(ii) = 1;
        nodePosNew(nNodesUnique,:) = pos;
        nodeDiamNew(nNodesUnique) = nodeDiam(ii);
        validatedNodesNew(nNodesUnique) = 1;
    end
    
end

% Close the wait bar.
close(hwait);

%%% Find the new edges based upon nodeMap
nodeEdges = nodeMap(nodeEdges);

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
% use the nodeMap to find when nB_new=1 and nB_old=2 for given nodes and
% then delete all nodes and edges back to the bifurcation node

%% Reassign output variables
nodes = nodePosNew;
validatedNodes = validatedNodesNew';
edges = nodeEdges;
fprintf('Regraph reduced %d nodes to %d\n',nNodes,size(nodes,1))

%% Assign the validated edges index
validatedEdges = zeros(size(edges,1),1);
for uu = 1:size(edges,1)
    if validatedNodes(edges(uu,1)) == 1 && validatedNodes(edges(uu,2)) == 1
        validatedEdges(uu) = 1;
    end
end

end