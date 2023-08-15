function [nodePos2,nodeDiam2,nodeDiamThetaIdx2,nodeBC2,nodeBCType2,nodeType2,nodeSegN2,nodeEdges2,edgeFlag2] = removeNodes( nodeFlag, nodePos, nodeDiam, nodeDiamThetaIdx, nodeBC, nodeBCType, nodeType, nodeSegN, nodeEdges );

hwait = waitbar(0,'Removing deleted nodes');

% remove abandoned nodes
nNodes = 0;
nodePosTmp = [];
nodeDiamTmp = [];
nodeDiamThetaIdxTmp = [];
nodeBCTmp = [];
nodeBCTypeTmp = [];
nodeTypeTmp = [];
nodeMap = zeros(size(nodePos,1),1);
for iN = 1:size(nodePos,1)
    if nodeFlag(iN)==1
        nNodes = nNodes + 1;
        nodePosTmp(nNodes,:) = nodePos(iN,:);
        nodeDiamTmp(nNodes) = nodeDiam(iN);
        nodeDiamThetaIdxTmp(nNodes) = nodeDiamThetaIdx(iN);
        nodeBCTmp(nNodes) = nodeBC(iN);
        nodeBCTypeTmp(nNodes) = nodeBC(iN);
        nodeTypeTmp(nNodes) = nodeType(iN);
        nodeSegNTmp(nNodes) = nodeSegN(iN);
        nodeMap(iN) = nNodes;
    end
end
nodePos2 = nodePosTmp;
nodeDiam2 = nodeDiamTmp;
nodeDiamThetaIdx2 = nodeDiamThetaIdxTmp;
nodeBC2 = nodeBCTmp;
nodeBCType2 = nodeBCTmp;
nodeType2 = nodeTypeTmp;
nodeSegN2 = nodeSegNTmp;
nodeEdges = nodeMap(nodeEdges);

% remove edges with abandoned nodes (i.e. = 0)
[ir,ic] = find(nodeEdges==0);
nodeEdges(ir,:) = [];

% remove redundant edges
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
[b,i,j]=unique(sE);
nodeEdges2 = nodeEdges(sort(i),:);

edgeFlag2 = zeros(1,size(nodeEdges,1));

close(hwait)
