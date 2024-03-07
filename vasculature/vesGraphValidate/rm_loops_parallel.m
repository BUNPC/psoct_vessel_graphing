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

%% Initialize graph
g = graph(edges(:,1), edges(:,2));
% Boolean for tracking whether loops exist
tf = hascycles(g);

%% While graph conatins loops, continue removing them in parallel
while tf
    % Separate into subgraphs, remove loops, recombine subgraphs
    [nodes, edges] =...
        subgraph_rm(nodes, edges, angio, delta, v_min, mv_iter, viz, g);

    % Determine whether updated graph contains cycles (hascycles)
    g = graph(edges(:,1), edges(:,2));
    tf = hascycles(g);
    
    % If tf is false, then there are no loops remaining in the graph, and
    % this function will return to parent call.
end


function [n, e] = subgraph_rm(nodes, edges, angio, delta, v_min, mv_iter, viz, g)
    %SUBGRAPH_RM create subgraphs, remove loops, recombine into graph
    
    % Debugging variables
    tmp = []; tmp2 = [];
    % Find loops in graph. Include a maximum to remain within memory limit
    [cnodes, cedges] = allcycles(g,'MaxNumCycles',500);
    cnodes_lst = horzcat(cnodes{:});

    %%% Find connected components (subgraphs)
    % This will classify nodes into the same "bin" index if they are
    % connected. The "weak" type ensures that connected branches all have
    % the same bin index.
    bins = conncomp(g);
   
    %%% Find subgraphs containing loops (compare cnodes and bins)
    % Create an array of subraph indices containing nodes in loops
    subgraph_idcs = bins(cnodes_lst);
    % Identify the unique subgraph indices
    sub_idcs_uq = unique(subgraph_idcs);
    
    %%% Separate graph into subgraphs
    
    
    %%% Store loop subgraphs in struct


    %%% Create parallel threads to remove loops
    
    
    %%% Recombine subgraphs (after loop removal) and subgraphs (w/o loops)

end

end