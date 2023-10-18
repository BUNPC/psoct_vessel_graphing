function [nodes, edges] = rm_loops(nodes, edges, angio, delta, v_min, mv_iter, lim)
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
%       delta (int): Search radius for x,y,z directions (units = voxels)
%       v_min (double): minimum voxel intensity threshold for the move to
%               mean function. The new voxel position will only be
%               reassigned if the voxel intensity of the new node position
%               is >= v_min.
%       mv_iter (int): number of iterations in move to mean.
%       lim (struct): structure with graph limits for visualization
%   OUTPUTS:
%       n ([n,3] array): node locations
%       e ([m,2] array): edges connecting each node

%% TODO
% 1) downsample_loops function increases the number of loops on some
% iterations when using noisy data.
%   - downsample each loop individually
%   - write logic to determine when a loop has been fully downsampled (ie,
%       number of nodes = number of edges)
%
% 2) Compare the function rm_reindex and my code to see if they are the same
%
% 3) Move to mean operates on all nodes of the graph, rather than just the
% nodes in the loops. Code will run faster if performed on just the loops.

%% Initialize variables
% Count number of loops in graph
[n_pre, cnodes, cedges] = count_loops(edges);
% Initialize maximum number of loops to first count of loops
n_max = n_pre;

% Flag to protect 
protect = true;

% Struct to store nodes, edges for move to mean function
im_mv = struct();
im_mv.angio = angio;

% Counter to track number of iterations
cnt = 1;

% Counter to track number of increments to delta (search radius)
dcnt = 0;

% Plotting limits for viewing cycles
lim.x = [70, 100];
lim.y = [60, 100];
lim.z = [60, 100];

% Graph prior to preprocessing
visualize_graph(nodes, edges, 'Before Loop Removal', []);

%% Determine if the loop is sparse
% Compare number of adjoined segments to number of nodes. If these
% are equal, then down sampling will not remove the node. The
% solution would be to remove the longest edge.

%%% Compute sparsity for each loop in the graph
% This returns a binary array (1=sparse, 0=not sparse)
loop_sparsity = graph_sparsity(edges);

% Indices of sparse loops
sidx = find(loop_sparsity); %#ok<EFIND> 

% If there are sparse loops
if ~isempty(sidx)
    % Subset of cell array of nodes in loop
    
    % Call function to remove longest edge from sparse loops
    edges_rm = rm_loop_edge(nodes, edges, nidx_loop);
    % Reassign output for continuity in function
    edges = edges_rm;
    % Recalculate loops
    [n_pre, cnodes, cedges] = count_loops(edges);
end

%% Perform move to mean, down sample, and remove longest edge
while ~isempty(cnodes)
    %%% Debugging information + Visualization
    % Outer while-loop iteration
    fprintf('\nLoop Removal Iteration = %i\n', cnt)    
    % Visualize Graph
    iter_str = strcat('Iteration ',num2str(cnt));
    tstr = {'Before Mv2Mean & Downsampling Graph',iter_str};
    visualize_graph(nodes, edges, tstr, []);
    xlim(lim.x); ylim(lim.y); zlim(lim.z);

    %% Initialize variables for inner while-loop
    % Nodes from first loop in list
    n_idcs = cnodes{1,:};
    % Edge indices from first loop in list
    e_idcs = cedges{1,:};
    % Number of loops before inner while-loop
    n_pre = length(cnodes);
    % Number of loops after down sample, move to mean, remove longest edge
    % Initialize this equal to n_pre, before entering inner while-loop
    npost = n_pre;
    while npost >= n_pre        
        %% Move to mean (collapse loops)
        % Create struct of graph to be compatible with move to mean
        im_mv.nodes = nodes;
        im_mv.edges = edges;
        % Perform move to mean to collapse nodes.
        for j=1:mv_iter
            im_mv = mv_to_mean(im_mv, v_min);
        end
        % Extract node positions + edges from struct (edges are unchanged)
        nodes_mv = im_mv.nodes;    
        edges_mv = im_mv.edges;
        % Visualize Graph
        tstr = {'After Mv2Mean',iter_str};
        visualize_graph(nodes_mv, edges_mv, tstr, []);
        xlim(lim.x); ylim(lim.y); zlim(lim.z);
    
        %% Regraph (downsample) to remove collapsed loops
        [nodes_ds, edges_ds] =...
            downsample_loops(n_idcs, nodes_mv, edges_mv, delta, protect);
        % Visualize after downsample
        tstr = {'After Mv2Mean + Downsample',iter_str};
        visualize_graph(nodes_ds, edges_ds, tstr, []);
        xlim(lim.x); ylim(lim.y); zlim(lim.z);

        %% Check for sparse loops. If exist, remove longest edge
        % TODO: this section replicates the one above. Place into fun.

        % Function to find sparse loops

        % Function to remove longest edge

        % Reassign to nodes_ds, edges_ds for continuity

        %% Determine if number of loops decreased
        % TODO: rethink logic here

        % Number of loops after downsampling
        [npost, cnodes, cedges] = count_loops(edges_ds);
        
        % If the number of loops reduced
        if npost < n_pre
            % Reassign nodes/edges and leave inner while-loop 
            edges = edges_ds;
            nodes = nodes_ds;
            break
        % If the loop is fully down sampled
        elseif 1==1
            edges = rm_loop_edge(nodes_ds, edges_ds);
            nodes = nodes_ds;
        else
            edges = edges_ds;
            nodes = nodes_ds;
        end

        % Recalculate number of loops post processing
        [npost, cnodes, cedges] = count_loops(edges);
    end
end

%% Initial loop-removal logic
%{
while ~isempty(cnodes)
    %%% Print iteration
    fprintf('\nLoop Removal Iteration = %i\n', cnt)
    
    % Visualize Graph
    iter_str = strcat('Iteration ',num2str(cnt));
    tstr = {'Before Mv2Mean & Downsampling Graph',iter_str};
    visualize_graph(nodes, edges, tstr, []);
    xlim(lim.x); ylim(lim.y); zlim(lim.z);

    %% Move to mean (collapse loops)
    % Create struct of graph to be compatible with move to mean
    im_mv.nodes = nodes;
    im_mv.edges = edges;
    % Perform move to mean to collapse nodes.
    for j=1:mv_iter
        im_mv = mv_to_mean(im_mv, v_min);
    end
    % Extract node positions + edges from struct (edges are unchanged)
    nodes_mv = im_mv.nodes;    
    edges_mv = im_mv.edges;
    % Visualize Graph
    tstr = {'After Mv2Mean',iter_str};
    visualize_graph(nodes_mv, edges_mv, tstr, []);
    xlim(lim.x); ylim(lim.y); zlim(lim.z);

    %% Regraph (downsample) to remove collapsed loops
    [nodes_ds, edges_ds] =...
        downsample_loops(cnodes, nodes_mv, edges_mv, delta, protect);
    % Visualize after downsample
    tstr = {'After Mv2Mean + Downsample',iter_str};
    visualize_graph(nodes_ds, edges_ds, tstr, []);
    xlim(lim.x); ylim(lim.y); zlim(lim.z);

    %%% Number of loops after down sampling
    [n_post, cnodes] = count_loops(edges_ds);
    % If no loops remain, then leave while loop
    if isempty(cnodes)
        break
    else
        fprintf('Number of loops after Mv2Mean & downsample = %d\n', n_post)
    end

    %% Determine if fully down-sampled loops remain
    % The function "downsample_loops" preserves the end nodes. This will
    % preserve Loops with three or more connecting segments, once the loops
    % have just one node per edge.  
    % TODO:
    % 1) add logic to determine when the loop has been fully downsampled
    % 2) should only remove the longest edge once the number of edges
    % equals the number of nodes in the loop. Once this condition is met,
    % all the edges will have been fully down sampled, so removing the
    % longest edge will open the loop. Otherwise, the loop may open along a
    % segment and then be closed on the next down sample iteration.

    % Set the value of maximum number of loops across this process
    n_max = max(n_max, n_post);
    
    %%% Check if (nloops == nedges) for all loops
    % Remove longest edge

    %%% If the number of loops remained constant
    if (n_post == n_pre) && (dcnt < 3)
        % Unprotect the end point nodes.
        protect = false;
        % Increase search radius to try to collapse all loops
        delta = delta + 50;
        % Counter for number of delta increases
        dcnt = dcnt + 1;
    %%% If the delta has incremented 3 times & loops remain
    elseif dcnt >= 3
        % Remove longest edge of loop
        edges = rm_loop_edge(nodes_ds, edges_ds);
        % Count the number of loops
        [n_post, cnodes] = count_loops(edges);
        % Display number of loops
        fprintf('Number of loops after removing longest edge = %d\n',...
            n_post)
    %%% If number of loops increased after down sampling
    elseif (n_post > n_pre) || (n_post > n_max)
        % TODO: determine better logic here.
        %{
        % Use the graph after mv2mean and before down sample.
        % Remove longest edge.
        edges = rm_loop_edge(nodes_mv, edges_mv);
        % Verify number of loops decreased
        [nloops_rm, ~] = count_loops(edges);
        assert( (nloops_rm < n_post) ,...
            'Removing longest edge did not reduce the number of loops');
        % If condition met, then reassign n_post and continue
        n_post = nloops_rm;
        %}
    end

    % Iterate the loop counter
    cnt = cnt + 1;

    % Update the number of loops before next iteration
    n_pre = n_post;

    % Close figures from iteration
%     close all;
end
%}

%%% Visualize graph after removing loops
visualize_graph(nodes, edges, 'After Loop Removal', []);

end

%% Function to count number of loops in graph
function [nloops, cnodes, cedges] = count_loops(edges)
%   INPUTS:
%       edges ([n,2] array): edges of graph
%   OUTPUTS:
%       nloops (int): number of loops
%       cnodes (cell array): each cell array contains node indices
%                       belonging to a loop

% Convert edges to graph
g = graph(edges(:,1), edges(:,2));

% Generate index of loops
[cnodes, cedges] = allcycles(g);

% Count number of loops
nloops = length(cnodes);

end









