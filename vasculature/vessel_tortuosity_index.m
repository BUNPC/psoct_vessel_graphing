function vi_dist=vessel_tortuosity_index(Graph,thresh)
%% find the tortuosity of vessel segments of a graph
% input:
%       G - graph struct data
%       thresh - threshold of vessel length, unit: micron
% output:
%       vi_dist - vector of vessel tortuosity index
% Author: Jiarui Yang
% 11/02/20
segInfo = Graph.segInfo;
seg_index = find(segInfo.segLen_um(:) >= thresh/10);
vi_dist = zeros(1,length(seg_index));
for i = 1:length(seg_index)
    l_graph = segInfo.segLen_um(seg_index(i));
    endnodes = segInfo.segEndNodes(seg_index(i),:);
    e1 = Graph.nodes(endnodes(1),:);
    e2 = Graph.nodes(endnodes(2),:);
    e_dist = sqrt((e1(1)-e2(1))^2+(e1(2)-e2(2))^2+(e1(3)-e2(3))^2);
    vi_dist(i)=l_graph/(e_dist*2);
end

end