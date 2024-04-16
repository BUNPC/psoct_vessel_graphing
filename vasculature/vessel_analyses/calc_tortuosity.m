function tort= calc_tortuosity(data)
% data is .mat variable

nodes = data.Graph.nodes;
endnodes = data.Graph.segInfo.segEndNodes;
len = data.Graph.segInfo.segLen_um;
nves = length(data.Graph.segInfo.segLen_um);
vox_dim = data.Graph.vox;

tort = zeros(nves, 1);
    for j=1:nves
        % convert segment end nodes to cartesian coordinate
        node1 = nodes(endnodes(j,1), :);
        node2 = nodes(endnodes(j,2), :);
        % Convert cartesian coordinate to offset in microns
        node1 = node1 .* vox_dim;
        node2 = node2 .* vox_dim;
        % Calcualte euclidean distance of segment
        euc = sqrt((node1(1) - node2(1)).^2 +...
                    (node1(2) - node2(2)).^2 +...
                    (node1(3) - node2(3)).^2);
        % Calculate tortuosity (single arc-chord ratio)
        tort(j) = len(j) ./ euc;
    end
    % Remove infinite tortuosity (loops)
    tort(tort==inf) = [];

end