[l,pos]=find(Data.Graph.segInfo.nodeSegN(:)==idx);
Graph=[];
Graph.nodes=Data.Graph.nodes(l,:);
Graph.edges=[];
for i=1:length(l)
    node_idx=l(i);
    edge_tmp=find(Data.Graph.edges(:,1)==node_idx);
    for j=1:length(edge_tmp)
        [edge_idx,~]=find(l(:)==Data.Graph.edges(edge_tmp(j),2));
        if ~isempty(edge_idx)
            Graph.edges=[Graph.edges;i edge_idx];
        end
    end
end