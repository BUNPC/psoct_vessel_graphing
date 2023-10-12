function [nodes, edges] = rm_loops(nodes, edges, angio, loop_flag, delta, v_min, mv_iter)
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
%       loop_flag (logical):
%               true = down sample loops with "downsample_loops"
%               false = down sample entire graph
%       delta (int): Search radius for x,y,z directions (units = voxels)
%       v_min (double): minimum voxel intensity threshold for the move to
%               mean function. The new voxel position will only be
%               reassigned if the voxel intensity of the new node position
%               is >= v_min.
%       mv_iter (int): number of iterations in move to mean.
%   OUTPUTS:
%       n ([n,3] array): node locations
%       e ([m,2] array): edges connecting each node

%% TODO
% Some loops remain, even after increasing the delta an unprotecting the
% end nodes. Need to investigate root cause.
%
% The function "move to mean" is currently operating on all nodes of the
% graph, rather than just the nodes in the loops. Prior tests found that
% the vascular morphology is retained when running move to mean on the
% entire graph. It is possible to modify move to mean to only operate on
% the loops, but this is a low priority.

%% Convert from nodes + edges into Matlab graph
% Create standard Matlab graph
g = graph(edges(:,1), edges(:,2));

% Detect loops in graph
[cnodes, ~] = allcycles(g);

% Count total number of nodes belonging to loops
n_pre = length(cnodes);

% Flag to protect 
protect = true;

%% While loops exist: regraph + move to mean
% Counter to track number of iterations
cnt = 1;
% Counter to track number of increments to delta (search radius)
dcnt = 0;

while ~isempty(cnodes)
    %%% Print iteration
    fprintf('\nLoop Removal Iteration = %i\n', cnt)
    
    %%% Regraph (downsample) to remove collapsed loops
    % Initialized validated nodes vector so that regraph code runs
    validated_nodes = zeros(size(nodes,1),1);
    % Call function to down sample
    if loop_flag
        [nodes, edges] =...
            downsample_loops(cnodes, nodes, edges, delta, protect);
    else
        [nodes, edges, ~, ~] =...
            regraphNodes_new(nodes, edges, validated_nodes, delta);
    end

    %%% Move to mean (collapse loops)
    % Create copy of regraphed graph to be compatible with move to mean
    im_mv = struct();
    im_mv.nodes = nodes;
    im_mv.edges = edges;
    im_mv.angio = angio;
    % Perform move to mean to collapse nodes.
    for j=1:mv_iter
        im_mv = mv_to_mean(im_mv, v_min);
    end
    % Extract nodes from struct (edges are constant)
    nodes = im_mv.nodes;
    
    %%% Detect loops in graph
    g = graph(edges(:,1), edges(:,2));
    [cnodes, ~] = allcycles(g);
    
    %%% Determine if fully down-sampled loops remain
    % The function "downsample_loops" preserves the end nodes, so loops are
    % preserved if they contain 3 or more end nodes. This code will
    % determine if the total number of nodes belonging to loops remains
    % constant after down sampling.

    % Number of loop nodes after last iteration
    n_post = length(cnodes);
    fprintf('Number of edgecycles (loops) after DS + M2M = %d\n', n_post)

    %%% Compare number of loop nodes before/after last iteration
    % If the number of nodes remains constant (and loops remain)
    if (n_pre == n_post) && (dcnt < 3)
        % Unprotect the end point nodes.
        protect = false;
        % Increase search radius to try to collapse all loops
        delta = delta + 50;
        % Counter for number of delta increases
        dcnt = dcnt + 1;
    % Elseif delta has been incremented to > 100
    elseif (n_pre == n_post) && (dcnt >= 3)
        % Remove longest edge of loop
        edges = rm_loop_edge(cnodes, nodes, edges);
        %%% Detect loops in graph
        g = graph(edges(:,1), edges(:,2));
        [cnodes, ~] = allcycles(g);
    else
        % Update n_pre for next iteration
        n_pre = n_post;
    end

    % Iterate the loop counter
    cnt = cnt + 1;

    % Close figures from iteration
    close all;
end

%%% Visualize graph after removing loops
visualize_graph(nodes, edges, 'After Loop Removal', []);

end