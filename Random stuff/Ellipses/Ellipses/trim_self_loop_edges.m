%% trim the edges
global Data
% Graph.nodes = nodes;
% Graph.edges = edges;

sameEdgeIdx = [];
for u = 1:size(Data.Graph.edges,1)
    if Data.Graph.edges(u,1) == Data.Graph.edges(u,2)
        sameEdgeIdx = [sameEdgeIdx; u];
    end
end
Data.Graph.edges(sameEdgeIdx,:) = [];

%% save 3D graph
global Data
img=Data.angioVolume;
img(img(:)~=254&img(:)~=252)=0;
img(img(:)==254|img(:)==252)=1;
MAT2TIFF(permute(img,[2 3 1]),'marching_Vis_trimmed.tif');

%% vis the histogram
binEdge=0:3:45;
figure;histogram(ves,binEdge);
hold on;
histogram(mar,binEdge);

%% trim the short segments
global Data
shortSegIdx=[];
for u = 1:length(Data.Graph.segInfo.segLen)
    if Data.Graph.segInfo.segLen(u)<20
        shortSegIdx = [shortSegIdx; u];
    end
end

nodeRmvIdx=[];
for u=1:length(Data.Graph.segInfo.nodeSegN)
    if ismember(Data.Graph.segInfo.nodeSegN(u),shortSegIdx) || Data.Graph.segInfo.nodeSegN(u)==0
        nodeRmvIdx=[nodeRmvIdx; u];
    end
end

edgeIdx=[];
for u=1:size(Data.Graph.edges,1)
    if ismember(Data.Graph.edges(u,1),nodeRmvIdx) && ismember(Data.Graph.edges(u,2),nodeRmvIdx)
        % if both nodes are not to be removed
        edgeIdx=[edgeIdx; u];
    end
end

% Data.Graph.nodes(nodeRmvIdx,:)=[];
Data.Graph.edges(edgeIdx,:)=[];