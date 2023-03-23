function vi_dist = vessel_tortuosity_index(Graph, thresh)
%% find the tortuosity of vessel segments of a graph
% input:
%       G - graph struct data
%       thresh - threshold of vessel length, unit: micron
% output:
%       vi_dist - vector of vessel tortuosity index
% Author: Jiarui Yang
% 11/02/20

%% variable initialization
% length of each segment (microns)
seglen_mat = Graph.segInfo.segLen_um;
% Matrix of end nodes of all segment
endnode_mat = Graph.segInfo.segEndNodes;

% Find segments with length greater than (1/10) * threshold
seg_index = find(Graph.segInfo.segLen_um(:) >= thresh./10);
% Initialize array for storing distance of each vessel
vi_dist = zeros(1, length(seg_index));

%% Calculate length / (2 * Euclidean Distance)
for i = 1:length(seg_index)
    disp(i)
    % length of current segment in micron
    seglen = seglen_mat(seg_index(i));
    % e = end node indices of current segment
    e = endnode_mat(seg_index(i),:);
    % Retrieve (x,y,z) coordinates for each node (e1, e2)
    e1 = Graph.nodes(e(1),:);
    e2 = Graph.nodes(e(2),:);
    
    % Calculate Euclidean distance of nodes
    e_dist = sqrt( (e1(1)-e2(1))^2 + (e1(2)-e2(2))^2 + (e1(3)-e2(3))^2);

    % Calculate tortuosity
    vi_dist(i) = seglen ./ (e_dist .* 2);
end

end