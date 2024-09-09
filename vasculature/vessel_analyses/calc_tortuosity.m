function tort = calc_tortuosity(data)
%CALC_TORTUOSITY calculate tortuosity of each segment
% This uses the tortuosity metric known as the arc:chord ratio (total
% length divided by Euclidean distance). The output is an array of
% tortuosity vectors.

% INPUTS
%   data (struct): data structure containing:
%       nodes: node positions [y,x,z] (units = voxels)
%       end nodes: end node positions [y,x,z] (units = voxels)
%       len (double array): length of each vessel (units = voxels)
%       nves (double): number of vessels
% OUTPUTS
%       tort (double vector): vector of tortuosity measurements for each
%           segment

nodes = data.Graph.nodes;
endnodes = data.Graph.segInfo.segEndNodes;
len = data.Graph.segInfo.segLen;
nves = length(data.Graph.segInfo.segLen);

% Initialize the tortuosity array
tort = zeros(nves, 1);
for j=1:nves
    % convert segment end nodes to cartesian coordinate
    node1 = nodes(endnodes(j,1), :);
    node2 = nodes(endnodes(j,2), :);
    % Calcualte euclidean distance of segment (units = voxels)
    euc = sqrt((node1(1) - node2(1)).^2 +...
                (node1(2) - node2(2)).^2 +...
                (node1(3) - node2(3)).^2);
    % Calculate tortuosity (single arc-chord ratio)
    tort(j) = len(j) ./ euc;
end
% Remove infinite tortuosity (loops)
tort(tort==inf) = [];

end