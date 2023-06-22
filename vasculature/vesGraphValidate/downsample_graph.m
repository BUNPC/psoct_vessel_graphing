function [nodes, edges, validatedNodes,validatedEdges] =...
    downsample_graph(nodes, edges, validatedNodes, vox_xy, vox_z)
%%% Down sample 
% INPUTS:
%       nodes (matrix): Nodes prior to downsample
%       edges (matrix): edges prior to downsample
%       validatedNodes (array): validated node indices
%       vox_xy (double): voxel dimensions (x,y)
%       vox_z (double):  voxel dimension (z)
% OUTPUTS:
%       nodes (matrix): 
%       edges (matrix): 
%       validatedNodes (array): validated node indices
%       validatedEdges (array): validated edges
%{
TODO:
- rewrite this function so that the outer for-loop iterates over the
segments. Within a segment, iterate over the nodes and compare the location
of nodes within the boundaries [hxy, hxy, hz].
%}



%% Initialization
% Set node position and edges
nodePos = nodes;
nodeEdges = edges;
nNodes = size(nodePos,1);
nodeDiam = zeros(nNodes,1);

% set voxel dimensions
hxy = vox_xy;
hz = vox_z;

%%% Create node maps
nNodesUnique = 1;
nodeMap = zeros(nNodes,1);
nodeUnique = zeros(nNodes,1);
% Set first value of arrays
nodeMap(1) = 1;
nodePosNew = nodePos(1,:);
nodeUnique(1) = 1;
validatedNodesNew(1) = validatedNodes(1);

for ii=2:nNodes
    % Set position equal to current node in list.
    pos = nodePos(ii,:);
    sprintf('Processed %i out of %i nodes',ii,nNodes)
    if validatedNodes(ii)==0
        lst = find(...
            pos(1)>=(nodePosNew(:,1)-hxy) & pos(1)<=(nodePosNew(:,1)+hxy) & ...
            pos(2)>=(nodePosNew(:,2)-hxy) & pos(2)<=(nodePosNew(:,2)+hxy) & ...
            pos(3)>=(nodePosNew(:,3)-hz)  & pos(3)<=(nodePosNew(:,3)+hz));
        if isempty(lst)
            nNodesUnique = nNodesUnique+1;
            nodeMap(ii) = nNodesUnique;
            nodeUnique(ii) = 1;
            nodePosNew(nNodesUnique,:) = pos;
            nodeDiamNew(nNodesUnique) = nodeDiam(ii);
            validatedNodesNew(nNodesUnique) = 0;
        else
            if length(lst)>1
                clear d
                for iLst=1:length(lst)
                    d(iLst) = norm(pos-nodePosNew(lst(iLst),:));
                end
                [~, closestNode] = min(d);
            else
                closestNode = 1;
            end
            nodeMap(ii) = lst(closestNode);
        end
    % nodeValidated(ii)==1
    else
        nNodesUnique = nNodesUnique+1;
        nodeMap(ii) = nNodesUnique;
        nodeUnique(ii) = 1;
        nodePosNew(nNodesUnique,:) = pos;
        nodeDiamNew(nNodesUnique) = nodeDiam(ii);
        validatedNodesNew(nNodesUnique) = 1;
    end    
end
% Remap edges and nodes
nodeEdgesNew = nodeMap(nodeEdges);
nodeEdges = nodeEdgesNew;

%% What does this section do???
% prune edges - still need to handle small loops
% point edges
nodeEdges = nodeEdges(find(nodeEdges(:,1)~=nodeEdges(:,2)),:);

% redundant edges
sE = cell(size(nodeEdges,1),1);

for ii=1:length(nodeEdges)
    if nodeEdges(ii,1)<nodeEdges(ii,2)
        sE{ii} = sprintf('%05d%05d',nodeEdges(ii,1),nodeEdges(ii,2));
    else
        sE{ii} = sprintf('%05d%05d',nodeEdges(ii,2),nodeEdges(ii,1));
    end
end

[~,i,~]=unique(sE);
nodeEdges = nodeEdges(sort(i),:);

% check for new dangling nodes (i.e. nodes with nB=1 that were nB=2)
% if we want to implement this, we just need to have nB_old and nB_new and
% use the nodeMap to find when nB_new=1 and nB_old=2 for given nodes and
% then delete all nodes and edges back to the bifurcation node

%%% Set outputs
nodes = nodePosNew;
edges = nodeEdges;

% Update validated nodes/edges
validatedNodes = validatedNodesNew';
validatedEdges = zeros(size(edges,1),1);
for uu = 1:size(edges,1)
    if validatedNodes(edges(uu,1)) == 1 && validatedNodes(edges(uu,2)) == 1
        validatedEdges(uu) = 1;
    end
end

fprintf('Regraph reduced %d nodes to %d\n',nNodes,size(nodes,1))

