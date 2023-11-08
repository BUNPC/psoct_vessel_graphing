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
% 1) In the outer while-loop, I select a subset of nodes to downsample.
% This subset is the first row in the cell array "cnodes". However, in the
% case that there are nested/connected loops, these nodes will be contained
% within multiple cell array rows. Therefore, the code is currently unable
% to process nested/connected loops because the variable "n_idcs" does not
% contain all the nodes within the nested structure. This is raising an
% error later in the code when trying to graph because nodes are being
% excluded.
% Solution:
%   - For each iteration, find the cell array rows with shared nodes.
%   - Combine these rows into a single array
%   - This also should be updated at the end of the inner while-loop:
%       "n_idcs = cnodes{1,:};"
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
cnt_outer = 1;

% Original delta size
delta0 = delta;

% Graph prior to preprocessing
% visualize_graph(nodes, edges, 'Before Loop Removal', []);

%% Remove longest edge of sparse loops
[n_pre, cnodes, cedges, edges] = open_sparse_loops(nodes, edges);

%% Perform move to mean, down sample, and remove longest edge
while ~isempty(cnodes)
    %%% Initialize variables for inner while-loop
    % TODO: find overlapping cell array nodes
    % Nodes from first loop in list
    n_idcs = cnodes{1,:};
    % Debugging line
    [n_idcs] = find_multiloop_nodes(cnodes);

    
    % Restore delta if it was incremented
    delta = delta0;
    % Edge indices from first loop in list
    e_idcs = cedges{1,:};
    % Number of loops before inner while-loop
    n_pre = length(cnodes);
    % npost = # loops after down sample, move to mean, remove longest edge.
    % Initialize equal to n_pre before entering inner while-loop.
    npost = n_pre;

    %%% Debugging information + Visualization
    % Initialize inner while-loop counter
    cnt_inner = 1;
    % Outer while-loop iteration
    fprintf('\nOuter while-loop iteration = %i\n', cnt_outer)    
    % Visualize Graph
    iter_str = strcat('Iteration ',num2str(cnt_outer));
    tstr = {'Before Mv2Mean & Downsampling Graph',iter_str};
    visualize_graph(nodes, edges, tstr, []);
    
    %%% Set graph limits based upon current cycle under investigation
    % Coordinates of all nodes
    lims = nodes(n_idcs,:);
    % Limits of each axis
    lim.x = [min(lims(:,1)) - 10, max(lims(:,1)) + 10];
    lim.y = [min(lims(:,2)) - 10, max(lims(:,2)) + 10];
    lim.z = [min(lims(:,3)) - 10, max(lims(:,3)) + 10];
    % Set limits in graphical display
    xlim(lim.x); ylim(lim.y); zlim(lim.z);


    while npost >= n_pre        
        %% Debugging information
        fprintf('\nInner while-loop iteration = %i', cnt_inner)
        fprintf('\nNumber of nodes before mv2mean & DS = %i', length(nodes))
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
        %%% Recalculate loops
        [~, cnodes, ~] = count_loops(edges_mv);
        nnodes = length(cnodes{1,:});
        fprintf('\nNumber of nodes in loop after mv2mean = %i', length(nodes_mv));
        %% Regraph (downsample) to remove collapsed loops
        [nodes_ds, edges_ds] =...
            downsample_loops(n_idcs, nodes_mv, edges_mv, delta, protect);
        % Visualize after downsample
        tstr = {'After Mv2Mean + Downsample',iter_str};
        visualize_graph(nodes_ds, edges_ds, tstr, []);
        xlim(lim.x); ylim(lim.y); zlim(lim.z);

        %% Check for sparse loops. If exist, remove longest edge
        [npost, cnodes, cedges, edges_ds] = open_sparse_loops(nodes_ds, edges_ds);
        
        % Reassign edges, nodes for next iteration of while-loop
        edges = edges_ds;
        nodes = nodes_ds;

        %% Determine if number of loops decreased
        
        %%% If the number of loops decreased
        if npost < n_pre
            % Verify that the updated list of loop nodes does not contain
            % the node indices from this current iteration of while loop.
            cnodes_all = horzcat(cnodes{:});
            assert(all(~ismember(n_idcs, cnodes_all)),...
                '\nThe function removed a loop but not the one under consideration.');
            
            % Iterate the outer while-loop counter
            cnt_outer = cnt_outer + 1;
            % Leave inner while-loop
            break
        else
            % Iterate the inner while-loop counter
            cnt_inner = cnt_inner + 1;
            % Increase delta
            delta = delta + 10;
            % Reassign loop node indices for next iteration. The indices of
            % the nodes in the loop may change during downsampling.
            % Therefore, this must be updated after each iteration.
            n_idcs = cnodes{1,:};
        end
    end
end

%% Initial loop-removal logic
%{
while ~isempty(cnodes)
    %%% Print iteration
    fprintf('\nLoop Removal Iteration = %i\n', cnt_outer)
    
    % Visualize Graph
    iter_str = strcat('Iteration ',num2str(cnt_outer));
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
    cnt_outer = cnt_outer + 1;

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
%count_loops return number of loops in graph, node indices, & edge indices
%   INPUTS:
%       edges ([n,2] array): edges of graph
%   OUTPUTS:
%       nloops (int): number of loops
%       cnodes (cell array): each cell contains node indices
%                       belonging to a loop
%       cedges (cell array): each cell contains the respective edge indices
%

% Convert edges to graph
g = graph(edges(:,1), edges(:,2));

% Generate index of loops
[cnodes, cedges] = allcycles(g);

% Count number of loops
nloops = length(cnodes);

end

%% Function to open a sparse loop edge
function [nloops, cnodes, cedges, edges] = open_sparse_loops(nodes, edges)
%open_sparse_loops: remove longest edge of sparse loop
%   INPUTS:
%       nodes ([n,3] array): nodes of graph
%       edges ([n,2] array): edges of graph
%   OUTPUTS:
%       nloops (int): number of loops
%       cnodes (cell array): each cell contains node indices
%                       belonging to a loop
%       cedges (cell array): each cell contains the respective edge indices

%%% Identify sparse loops
sp = graph_sparsity(edges);
% Convert sparsity array to boolean
sp = boolean(sp);
% Count number of loops before removing
[npre, ~, ~] = count_loops(edges);

%%% Remove longest edge of sparse loops
if any(sp)
    % Generate graph
    g = graph(edges(:,1), edges(:,2));
    % Find cycles in graph
    [cnodes, cedges] = allcycles(g);
    % Keep node indices from sparse cycles
    cnodes(~sp) = [];
    cedges(~sp) = [];
    % Remove the longest edge from each sparse loop
    edges = rm_loop_edge(nodes, edges, sp, cnodes, cedges);
end

%%% Recalculate loops
[nloops, cnodes, cedges] = count_loops(edges);
%%% Debugging information
nrm = nloops - npre;
fprintf('\nSparse Loops Removed Before Down Sampling = %i\n', nrm)

end

%% Function to return all nodes in nested loop structure
function [idcs] = find_multiloop_nodes(cnodes)
%find_multiloop_nodes
% The cell array "cnodes" contains a row containing the node indices for
% a graph cycle. However, in the case that there are nested/connected
% loops, these nodes will be contained within multiple cell array rows.
% This function takes the node indices in the first row and then searches
% for rows containing the same node indices.
%
% In the event no other rows contains the node indices, it will only return
% the node indices in the first row. If another row contains at least one
% of the indices, then it will return the corresponding row(s). This will
% be repeated for the new row until all rows of the multiloop structure are
% returned.
%
% INPUTS:
%       cnodes (cell array): each cell contains node indices belonging to
%                               a loop
% OUTPUTS:
%   n_idcs (double array):

%%% Take the first entry of the cell array
idcs = cnodes{1,:};

%%% Initialize variable for comparing 
cellidx_past = zeros(length(cnodes),1);

% Find cell arrays with the same node indices
while 1
    %%% Find all cell arrays containing nodes from the first row (idcs)
    % This will output an [X, Y] logical matrix for each cell array row,
    % where X is the length of idcs and Y is the number of elements in row.
    % This logical matrix must be reduced twice to a single logical (0,1),
    % where 0 indicates no matching node index and 1 indicates at least 1
    % matching node index.
    cellidx = cellfun(@(x) x==idcs(:), cnodes, 'UniformOutput', false);
    % Reduce dimensionality (2D -> 1D)
    cellidx = cellfun(@(x) any(x), cellidx, 'UniformOutput', false);
    % Reduce dimensionality (1D -> single logical)
    cellidx = cell2mat(cellfun(@(x) any(x), cellidx, 'UniformOutput', false));

    %%% Determine if at least one matching node index
    % If at least one cell array row contains a matching index
    if any(cellidx) && any(cellidx_past ~= cellidx)
        % Retrieve cell array nodes from logical 1 rows
        extra = cnodes(cellidx);
        
        % Convert from cell array to 1D double array
        e = [];
        for ii = 1:length(extra)
            e = [e, cell2mat(extra(ii,:))];
        end
        
        % Add all unique node indices to idcs
        idcs = unique([idcs, unique(e)]); %#ok<AGROW> 

        % Store the logical array for comparison in next iteration
        cellidx_past = cellidx;

    % No additional cell array rows contain the nodes.
    else      
        % Exit the while loop and return "idcs" to outer function
        break
    end

end

end






















