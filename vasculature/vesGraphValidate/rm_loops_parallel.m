function [nodes, edges] = rm_loops_parallel(nodes, edges, angio, delta, v_min, mv_iter, viz)
%RM_LOOPS_PARALLEL Separate graph into subgraphs, create parallel thread
%for each subgraph, and remove loops from subgraph in each thread.
%   The "rm_loops" function removes loops in series. This is time consuming
%   when a graph contains a large amount of loops. This function will
%   call "allcycles" to identify the loops in a graph. Then, it will
%   subdivide the graph into disjoint subgraphs of connected components.
%   Finally, it will use parallel threads to remove loops from each of
%   these subgraphs.
%   INPUTS:
%       nodes ([n,3] array): node locations
%       edges ([m,2] array): edges connecting each node
%       angio (double matrix): PS-OCT intensity volume (vessels are
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
%
%{
% TODO:
- determine order of operations for removing nodes in loops:
    option 1:
        - find loop nodes in entire graph:
            [cnodes, cedges] = allcycles(g,'MaxNumCycles',500);
        - run remove_reindex_nodes with both cnodes, nodes, and edges
    option 2:
        - separate graph into subgraphs
        - check each subgraph for loops "hascycles"
        - flag these to remove loops
%}

%% Initialize graph
g = graph(edges(:,1), edges(:,2));
% Boolean for tracking whether loops exist
tf = hascycles(g);

%% While graph conatins loops, continue removing them in parallel
while tf
    % Separate into subgraphs, remove loops, recombine subgraphs
    [nodes, edges] = ...
        subgraph_rm(nodes, edges, angio, delta, v_min, mv_iter, viz, g);

    % Determine whether updated graph contains cycles (hascycles)
    g = graph(edges(:,1), edges(:,2));
    tf = hascycles(g);
    
    % If tf is false, then there are no loops remaining in the graph, and
    % this function will return to parent call.
end

sprintf('Finished Removing Loops')


function [nodes, edges] =...
        subgraph_rm(nodes, edges, angio, delta, v_min, mv_iter, viz, g)
    %SUBGRAPH_RM create subgraphs, remove loops, recombine into graph
    %%% Find connected components (subgraphs)
    % This will classify nodes into the same "bin" index if they are
    % connected in the graph. The "weak" type ensures that connected
    % branches all have the same bin index.
    bins = conncomp(g);
    % The unique number of bins is the total number of subgraphs
    subgraph_idcs = unique(bins);
    
    %%% Separate graph into subgraphs (re-index nodes/edges)
    subgraphs = struct();
    % Create subgraph for each bin index
    parfor ii = 1:length(subgraph_idcs)
        % Find node indices with bin value of ii. In "nodeidx" the nodes
        % with a bin value equal to ii will be set to 1 and the others will
        % be set to 0.
        nodeidx = bins == ii;

        % Create list of node indices to remove. These have a bin value
        % different from ii.
        rm_idcs = 1:size(nodes,1);
        rm_idcs = rm_idcs(~nodeidx);
        
        % Reindex nodes/edges for subgraph
        [nodes_sub, edges_sub] = remove_reindex_nodes(rm_idcs, nodes, edges);
        n_nodes = size(nodes_sub,1);
        n_edges = size(edges_sub,1);

        % Check if subgraph contains loops
        subg = graph(edges_sub(:,1), edges_sub(:,2));
        loop_tf = hascycles(subg);

        % Place nodes/edges into struct
        subgraphs(ii).nodes = nodes_sub;
        subgraphs(ii).edges = edges_sub;
        subgraphs(ii).loop_tf = loop_tf;
        subgraphs(ii).n_nodes = n_nodes;
        subgraphs(ii).n_edges = n_edges;
    end
    
    
    %%% identify subgraph w/ loops (from struct)
    sub_idx = [subgraphs.loop_tf].';
    % struct index of subgraph containing loops
    sub_idx = find(sub_idx);
    L = length(sub_idx);
    
    %% Remove loops from subgraphs
    % Initialize new struct for storing subgraphs after loop removal
    subg_rm = subgraphs(sub_idx);

    parfor j = 1:L
        % Extract subgraph nodes/edges
        sub_nodes = subg_rm(j).nodes;
        sub_edges = subg_rm(j).edges;

        % Remove loops in subgraph
        [nodes_sub_rm, edges_sub_rm] = ...
            rm_loops(sub_nodes, sub_edges, angio, delta, v_min, mv_iter, viz);
        
        % Find number of nodes and edges
        n_nodes = size(nodes_sub_rm, 1);
        n_edges = size(edges_sub_rm, 1);
        
        % Verify loops were removed
        subg = graph(edges_sub_rm(:,1), edges_sub_rm(:,2));
        loop_tf = hascycles(subg);
        
        % Add to struct
        subg_rm(j).nodes = nodes_sub_rm;
        subg_rm(j).edges = edges_sub_rm;
        subg_rm(j).loop_tf = loop_tf;
        subg_rm(j).n_nodes = n_nodes;
        subg_rm(j).n_edges = n_edges;

        % Print to console that loop was removed
        fprintf('\n\nSubgraph %i removed\n', j);
    end
    
    %% Recombine subgraphs (after loop removal) and subgraphs (w/o loops)
    
    %%% Move subg_rm back into subgraphs struct
    for ii = 1:L
        subgraphs(sub_idx(ii)) = subg_rm(ii);
    end
    
    %%% Initialize arrays for combining all nodes/edges from all subgraphs
    % into a single graph
    n_nodes = sum([subgraphs.n_nodes]);
    n_edges = sum([subgraphs.n_edges]);
    nodes = zeros(n_nodes, 3);
    edges = zeros(n_edges, 2);

    % Init variables for calculating the offset for node_rm
    n_idx0 = 1;
    e_idx0 = 1;
    e_offset = 0;

    % Iterate over subgraphs
    for ii = 1:length(subgraph_idcs)
        % Extract nodes/edges from subgraph ii
        nodes_sub = subgraphs(ii).nodes;
        edges_sub = subgraphs(ii).edges;

        % Check if the subgraph was completely removed (due to being a self
        % loop without any non-loop nodes)
        if isempty(nodes_sub)
            continue
        end
        
        % Calculate end node indices
        n_idxf = n_idx0 + size(nodes_sub, 1) - 1;
        e_idxf = e_idx0 + size(edges_sub, 1) - 1;
        
        % Add subgraph nodes and edges to main list for entire graph
        nodes(n_idx0:n_idxf,:) = nodes_sub;
        edges(e_idx0:e_idxf,:) = edges_sub + e_offset;
        
        % Increment array indices
        n_idx0 = n_idxf + 1;
        e_idx0 = e_idxf + 1;
        
        % Increment edge offset by the number of nodes in the current
        % subgraph.
        e_offset = e_offset + size(nodes_sub,1);
    end

    % Visualize graph
    if viz
        visualize_graph(nodes, edges, 'Graph after loop removal',[])
    end
end

end



















