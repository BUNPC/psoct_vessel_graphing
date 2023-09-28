function [nodes, edges] = rm_loops(nodes, edges, angio, v_min, loop_flag, delta)
%rm_loops Remove loops in graph.
%   Outline:
%       - Use graph function "allcycles" to find loops
%       - While loops exist:
%           - convert from graph back to nodes/edges
%           - regraph (downsample nodes within search radius)
%           - move to mean (collapse loops)
%
%   INPUTS:
%       nodes ([n,3] array): node locations
%       edges ([m,2] array): edges connecting each node
%       im.angio (double matrix): PS-OCT intensity volume (vessels are
%               bright)
%       v_min (double): minimum voxel intensity threshold. The new voxel
%               position will only be reassigned if the voxel intensity of
%               the new node position is >= v_min.
%       loop_flag (logical):
%               true = down sample loops with "downsample_subset"
%               false = down sample entire graph
%       delta (int): Search radius for x,y,z directions (units = voxels)
%   OUTPUTS:
%       n ([n,3] array): node locations
%       e ([m,2] array): edges connecting each node

%% TODO
% The function "downsample_subset" preserves the end nodes, so
% loops are preserved if they contain 3 or more end nodes.
% - Add an if-statement to "rm_loops" to check if this condition occurs.
%       Track the number of nodes in cycles for each iteration. Once this
%       number reaches a constant, then:
%           - set "rm_end_node = True"
%           - pass "rm_end_node" into "downsample_subset"
% - Update "downsample_subset" to include/exclude end nodes from the
%       "pos_new" array of unique nodes.

%% Convert from nodes + edges into Matlab graph
% Create standard Matlab graph
g = graph(edges(:,1), edges(:,2));

% Detect loops in graph
[cnodes, ~] = allcycles(g);

%% While loops exist: regraph + move to mean
while ~isempty(cnodes)
    %%% Regraph (downsample) to remove collapsed loops
    % Initialized validated nodes vector so that regraph code runs
    validated_nodes = zeros(size(nodes,1),1);
    % Call function to regraph
    if loop_flag
        [nodes, edges] = downsample_subset(cnodes, nodes, edges, delta);
    else
        [nodes, edges, ~, ~] =...
            regraphNodes_new(nodes, edges, validated_nodes, delta);
    end

    %%% Move to mean (collapse loops)
    % Create copy of regraphed graph
    im_mv = struct();
    im_mv.nodes = nodes;
    im_mv.edges = edges;
    im_mv.angio = angio;
    % Perform move to mean five times to collapse nodes. This number was
    % determined emprically from testing a data subset.
    for j=1:5
        im_mv = mv_to_mean(im_mv, v_min);
    end
    
    %%% Detect loops in graph
    g = graph(edges(:,1), edges(:,2));
    [cnodes, ~] = allcycles(g);
end

end