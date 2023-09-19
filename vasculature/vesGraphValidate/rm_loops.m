function [nodes, edges] = rm_loops(nodes, edges, angio, v_min)
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
%   OUTPUTS:
%       n ([n,3] array): node locations
%       e ([m,2] array): edges connecting each node

%% Convert from nodes + edges into Matlab graph
% Create standard Matlab graph
g = graph(edges(:,1), edges(:,2));

% Detect loops in graph
cycles = allcycles(g);

%% While loops exist: regraph + move to mean
while ~isempty(cycles)
    %%% Regraph (downsample) to remove collapsed loops
    % Initialized validated nodes vector so that regraph code runs
    validated_nodes = zeros(size(nodes,1),1);
    % Search delta for the x,y,z directions (units = voxels)
    delta = 2;
    % Call function to regraph
    [nodes, edges, ~, ~] =...
        regraphNodes_new(nodes, edges, validated_nodes, delta);
    
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
    
    %%% Regraph to smooth collapse loops
    [nodes, edges, ~, ~] =...
        regraphNodes_new(im_mv.nodes, im_mv.edges, validated_nodes, delta);
    
    %%% Detect loops in graph
    g = graph(edges(:,1), edges(:,2));
    cycles = allcycles(g);
end

end