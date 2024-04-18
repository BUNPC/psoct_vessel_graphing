function [nodes, edges] = rm_loops(nodes, edges, angio, delta, v_min, mv_iter, viz)
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
%       viz (boolean): true = display figures of graph
%   OUTPUTS:
%       n ([n,3] array): node locations
%       e ([m,2] array): edges connecting each node

%% TODO
% 2) Compare the function rm_reindex and my code to see if they are the same

%% Initialize variables
% Flag to protect the nodes connected to the loops
protect = true;

% Struct to store nodes, edges for move to mean function
im_mv = struct();
im_mv.angio = angio;

% Counter to track number of iterations
cnt_outer = 0;

% Original delta size
delta0 = delta;

% Graph prior to preprocessing
if viz
    visualize_graph(nodes, edges, 'Before Loop Removal', []);
end

%% Remove longest edge of sparse loops
[~, cnodes, ~, edges] = open_sparse_loops(nodes, edges, viz);

%% Perform move to mean, down sample, and remove longest edge
while ~isempty(cnodes)
    %%% Initialize variables for inner while-loop
    % Find indices of first loop and any connected loop
    [n_idcs, n_loops] = find_multiloop_nodes(cnodes);    
    % Reset delta to original increment
    delta = delta0;
    % Number of loops in graph
    n_pre = length(cnodes);
    % Initialize number of loops after loop removal
    n_post = n_pre;

    %%% Debugging information + Visualization
    close all;
    % Iterate outer-while loop
    cnt_outer = cnt_outer + 1;
    % Initialize inner while-loop counter
    cnt_inner = 1;
    % Visualize Graph
    iter_str = strcat('Iteration ',num2str(cnt_outer));
    tstr = {'Before Mv2Mean & Downsampling Graph',iter_str};
    if viz
        visualize_graph(nodes, edges, tstr, []);
    end
    
    %%% Set graph limits based upon current cycle under investigation
    % Coordinates of all nodes
    lims = nodes(n_idcs,:);
    % Limits of each axis
    lim.x = [min(lims(:,1)) - 10, max(lims(:,1)) + 10];
    lim.y = [min(lims(:,2)) - 10, max(lims(:,2)) + 10];
    lim.z = [min(lims(:,3)) - 10, max(lims(:,3)) + 10];
    % Set limits in graphical display
    xlim(lim.x); ylim(lim.y); zlim(lim.z);

    while (n_pre - n_loops) < n_post
        %% Move to mean (collapse loops)
        % Create struct of graph to be compatible with move to mean
        im_mv.nodes = nodes;
        im_mv.edges = edges;
        for j=1:mv_iter
            % Move the loop nodes to the mean
            im_mv = mv_to_mean(im_mv, v_min, n_idcs);
        end
        % Extract node positions + edges from struct (edges are unchanged)
        nodes_mv = im_mv.nodes;    
        edges_mv = im_mv.edges;
        % Visualize Graph
        tstr = {'After Mv2Mean',iter_str};
        if viz
            visualize_graph(nodes_mv, edges_mv, tstr, []);
        end
        xlim(lim.x); ylim(lim.y); zlim(lim.z);
        %% Regraph (downsample) to remove collapsed loops
        [nodes_ds, edges_ds] =...
            downsample_loops(n_idcs, nodes_mv, edges_mv, delta, protect);
        % Visualize after downsample
        tstr = {'After Mv2Mean + Downsample',iter_str};
        if viz
            visualize_graph(nodes_ds, edges_ds, tstr, []);
        end
        xlim(lim.x); ylim(lim.y); zlim(lim.z);

        %% Check for sparse loops. If exist, remove longest edge
        % n_post is the number of loops after the inner while loop
        [n_post, cnodes, ~, edges_ds] = open_sparse_loops(nodes_ds, edges_ds, viz);
        
        % Reassign edges, nodes for next iteration of while-loop
        edges = edges_ds;
        nodes = nodes_ds;

        %% Update counters for next iteration
        % Iterate the inner while-loop counter
        cnt_inner = cnt_inner + 1;
        % Increase delta
        delta = delta + 4;
        % Reassign loop node indices for next iteration. The indices of
        % the nodes in the loop may change during downsampling.
        % Therefore, this must be updated after each iteration.
        if ~isempty(cnodes)
            n_idcs = cnodes{1,:};
        end
    end
    % Graph After Loop Removed
    tstr = {'Loop Removed',iter_str};
    if viz
        visualize_graph(nodes, edges, tstr, []);
    end
end

%%% Visualize graph after removing loops
if viz
    visualize_graph(nodes, edges, 'After Loop Removal', []);
end
%%% Close all figures prior to next iteration
close all;

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
[cnodes, cedges] = allcycles(g,"MaxNumCycles",500);

% Count number of loops
nloops = length(cnodes);

end

%% Function to open a sparse loop edge
function [nloops, cnodes, cedges, edges] = open_sparse_loops(nodes,edges,viz)
%open_sparse_loops: remove longest edge of sparse loop
%   INPUTS:
%       nodes ([n,3] array): nodes of graph
%       edges ([n,2] array): edges of graph
%       viz (bool): true = display debugging graph figures
%
%   OUTPUTS:
%       nloops (int): number of loops
%       cnodes (cell array): each cell contains node indices
%                       belonging to a loop
%       cedges (cell array): each cell contains the respective edge indices

%%% Identify sparse loops
[sp, cnodes] = graph_sparsity(edges);
% Convert sparsity array to boolean
sp = logical(sp);

%%% Remove longest edge of sparse loops
if any(sp)
    % Keep node indices from sparse cycles
    cnodes(~sp) = [];
    % Remove the longest edge from each sparse loop
    edges = rm_loop_edge(nodes, edges, sp, cnodes, viz);
end

%%% Recalculate loops
[nloops, cnodes, cedges] = count_loops(edges);

end

%% Function to return all nodes in nested loop structure
function [idcs, n_loops] = find_multiloop_nodes(cnodes)
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
cellidx_past(1) = 1;
combine_flag = 1;

%%% Find cell arrays with the same node indices
while combine_flag
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
        % Calculate total number of loops
        n_loops = sum(cellidx(:));
        % Exit the while loop and return "idcs" to outer function
        combine_flag = 0;
    end
end

end