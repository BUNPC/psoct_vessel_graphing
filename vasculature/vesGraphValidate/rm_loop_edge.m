function [edges_rm] = rm_loop_edge(nodes, edges, sp, nsp, viz)
%rm_loop_edge Remove longest edge of loop.
% -------------------------------------------------------------------------
% PURPOSE:
% This function handles the case where there are more than two segments
% involved in a loop. In this case, downsampling will just replace each
% endpoint with another node further away from the loop, thus expanding the
% loop. In this case, the longest edge will be removed, removing the loop.
% -------------------------------------------------------------------------
% INPUTS
%   g (graph): the original graph struct
%   nodes (array): [X, Y, Z] matrix of coordinates for nodes
%   edges (array): [M, 2] matrix of edges connecting node indices
%   nsp (cell array): Sparse loop nodes. Each row contains node indices in
%                       a sparse loop.
%	viz (bool): 1 = visualize graph at intermediate steps
%
% OUTPUTS
%   edges_rm (matrix): edge matrix without longest edge for each loop
%
% ISSUE: the matlab function "allcycles" is supposed to return edge indices
% belonging to cycles in a graph. However, it consistently returns the
% incorrect edge indices. There is no documentation online regarding this
% issue.
% RESOLUTION:
%   - write custom function to determine edges belonging to loops
%   - Incorporate into this function

%% Remove longest edge from one cycle until all sparse cycles are opened
cnt = 1;
while any(sp)

    %% Extract the first sparse loop nodes and edges from cell array
    % Node indices
    cnodes = nsp{1,:};
    
    % Find the edge indices where both start/end node belong to loops
    tf = ismember(edges, cnodes);
    tf = tf(:,1) & tf(:,2);
    
    % Extract the loop edges matrix [source node, target node]
    cedges = edges(tf,:);

    %%% Set graph limits based upon current cycle under investigation
    % Coordinates of all nodes
    lims = nodes(cnodes,:);
    % Limits of each axis
    lim.x = [min(lims(:,1)) - 10, max(lims(:,1)) + 10];
    lim.y = [min(lims(:,2)) - 10, max(lims(:,2)) + 10];
    lim.z = [min(lims(:,3)) - 10, max(lims(:,3)) + 10];

    %% Calculate distances for each edge
    % Initialize matrix for storing distances
    dmat = zeros(1,size(cedges,1));
    
    % Iterate over each edge and calculate distance b/w nodes
    for n = 1:length(dmat)
        % Extract nodes for edge
        e = cedges(n,:);
        n1 = nodes(e(1),:);
        n2 = nodes(e(2),:);
        % Calculate distance
        dmat(n) = pdist2(n1, n2);
    end
    
    %% Remove long edge from "edges" matrix

    %%% Find end nodes of longest edge
    % Find dmat index of longest edge (longest eucl. dist. to other nodes)
    [~, longest_idx] = max(dmat,[],"all","linear");
    % Convert from dmat index (same as cedges index) to start/end nodes
    e = cedges(longest_idx,:);
    
    %%% Remove longest edge from "edges" matrix
    % Find edge containing start/end nodes
    e_idx = (edges(:,1)==e(1) & edges(:,2)==e(2));
    % Remove edge from edges matrix
    edges(e_idx,:) = [];

    %% Recalculate graph sparsity for while-loop
    sp = graph_sparsity(edges);
    % Generate graph
    g = graph(edges(:,1), edges(:,2));
    % Find cycles in graph
    [nsp, ~] = allcycles(g);
    % Convert sparsity array to boolean
    sp = logical(sp);
    % Keep node indices from sparse cycles
    nsp(~sp) = [];
    
    %% Plot updated graph    
    tstr = strcat('Iteration ', num2str(cnt));
	if viz
		visualize_graph(nodes, edges, {'Longest Edge Removed',tstr},[]);
	end
    % Set limits in graphical display
    xlim(lim.x); ylim(lim.y); zlim(lim.z);
    cnt = cnt+1;
end


%%% Reassign output variable to make Matlab happy
edges_rm = edges;

end

