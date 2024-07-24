function seg_diam = calc_segment_diameter(node_seg, angio, nodes, ithresh, vox_dim)
% Calculate median diameter for each segment. 
%   Extract the nodes for each segment. Measure the diameter at each node
%   within each segment. Retain the median diameter for each segment. This
%   function calls the function "calc_diameter," which is part of the
%   imview3D suite.
%
%   INPUTS:
%       node_seg (uint8 matrix) - segment index for each node. Each node
%           has a unique segment index. This is the substruct:
%           "Data.Graph.segInfo.nodeSegN," where the "Data" struct is the
%           output of the function "init_graph."
%       angio (logical matrix) - PSOCT volume or segmentation volume
%       nodes (double matrix) - node position from graph (in voxels)
%       edges (double matrix) - edge indices from graph
%       ithresh (double) - approximate threshold of the segmentation
%       vox_dim (double [1,3]) - voxel size (microns)
%          
%   OUTPUTS:
%       seg_diam (array) - Median diameter of each segment (microns).

%% Initialize variables
% Identify all unique semgent indices
unique_segs = unique(node_seg);
% Initialize array to store median diameters
seg_diam = zeros(length(unique_segs),1);

%% Iterate over segments
for ii = 1:length(unique_segs)
    % Retrieve ii segment index
    seg_n = unique_segs(ii);
    % Extract nodes from segment index
    idx = node_seg == seg_n;
    seg_nodes = nodes(idx,:);
    % Function to calculate diameter for all nodes
    d = calc_diameter(angio, seg_nodes, ithresh, vox_dim);
    % Store the median diameter
    seg_diam(ii) = median(d);
end
