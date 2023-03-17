%% Prep environment & load the data w/ pruned segments
% Data.Graph.edges is a [N, 2] matrix, where each row contains the
% indices that the edge connects.
% Data.Graph.nodes is a [N, 3] matrix, where each row contains the [x,y,z]
% coordinates of the node.
% Data.Graph.nB is a 1D array containing the number of bifrucations at each
% end node.


clear; close all; clc;
tmp = load('Data_segpruned.mat');
Data = tmp.Data;
% Top-level directory
prel = strcat(pwd, '\GrTheory\');
addpath prel

% Create local variable of number of bifurcations for each end node
nb = Data.Graph.nB;
% Create local variable of all edges in Graph
edges = Data.Graph.edges;

%% Identify 
% Iterate through all branching nodes, identify the loops, & delete 1 edge
% if it meets the condition.

nnodes = size(Data.Graph.nodes,1);
lstN = 1:nnodes;

% Create list of nodes with > 2 bifurcations
bif3 = find(nb > 2); % nB = number of bifurcations
nbif3 = length(bif3);
hwait = waitbar(0,'Pruning loops');
n1 = 0; % number of nodes for which no loops removed
n2 = 0; % number of loops not pruned for skipped nodes
while ~isempty(bif3)
    bif3_end = bif3(end);
    bif3 = bif3(1:end-1);
    waitbar((nbif3-length(bif3))/nbif3, hwait);

    % [x,y,z] position of node w/ 3 bifurcations
    bif3_pos = Data.Graph.nodes(bif3_end,:);

    % Calculate Euclidean distances b/w each node and last node in bif3 array
    d = sum([...
        (Data.Graph.nodes(:,1) - ones(nnodes,1)*bif3_pos(1)) ...
        (Data.Graph.nodes(:,2)-ones(nnodes,1)*bif3_pos(2))...
        (Data.Graph.nodes(:,3)-ones(nnodes,1)*bif3_pos(3)) ].^2, 2 ).^0.5;
    
    %% Find which Euclidean distances < 30 (what units?, why use 30?)
    % Each iteration of the while loop selects the last element of the
    % array bif3. The Euclidean distance is calculated between this bif3
    % node and all other nodes. If the Euclidean distance is < 30, then the
    % index of the node's index will be saved to the array prox_node_idx .
    dthresh = 30;
    [prox_node_idx ,~] = find(d < dthresh);
    
    %% Find Graph.edges connected to bif3
    % Iterate through the list (prox_node_idx) of indices of nodes < 30 units
    % from bif3_end.
    %
    % Find the corresponding indices in edges. These indices in
    % edges correspond to the edges connecting to at least one
    % of these indices.
    % 
    % edges is a [N, 2] matrix, where each row contains the
    % indices that the edge connects.
    prox_edge_idx = [];
    for jj = 1:length(prox_node_idx)
        [ir, ~] = find(edges == prox_node_idx(jj));
        prox_edge_idx = [prox_edge_idx; ir];
    end
    % Remove redundant indices
    prox_edge_idx = unique(prox_edge_idx);
    
    %% Copy & re-index prox_edge_idx
    % The array prox_node_idx contains the absolute nodes indices from
    % Data.Graph.nodes. This section creates a local variable "prox_edges" which
    % contains the edges within proximity of the bif3 node. These edges are
    % indexed to the absolute index from Data.Graph.nodes.
    %
    % This loop will re-index prox_edges to the index in the array prox_node_idx

    prox_edges = edges(prox_edge_idx,:);
    for jj = 1:length(prox_node_idx)
        prox_edges( prox_edges == prox_node_idx(jj) ) = jj;
    end
    
    %% Find dangling nodes in "prox_edges"
    
    % Locate the proximal edges with vertices outside the proximity radius
    % of bif3. Confused as to how this loop handles this situation.
    % Why does it compare proximal edges to the length of the array rather
    % than the to length of the proximity radius?
    lst = find(prox_edges > length(prox_node_idx));
    for jj = 1:length(lst)
        kk = find( prox_edges(lst(jj)) == prox_node_idx );
        if isempty(kk)
            prox_node_idx(end+1) = prox_edges(lst(jj));
            prox_edges(lst(jj)) = length(prox_node_idx);
        else
            prox_edges(lst(jj)) = kk;
        end
    end
    %}
    %% Find edges connected to nodes with 3+ edges
    % edge3plus = edges connecting two nodes satisfying these conditions:
    %       a) contain > 2 bifurcations
    %       b) within 30 units of bif3
    
    % Find indices within range of endnodes (nB)
    ref1 = prox_node_idx(prox_edges(:,1));
    ref2 = prox_node_idx(prox_edges(:,2));
    % Remove indices exceeding the bounds of endnodes
    ref1( ref1 > length(nb) ) = [];
    ref2( ref2 > length(nb) ) = [];
    % Extract end nodes with >2 bifurcations & within 30 units of bif3
    edge_nb1 = nb( ref1 );
    edge_nb2 = nb( ref2 );
    
    %%% Find which edges connect to 2 nodes satisfying all above conditions
    % Pad shorter array with zeros to allow for boolean logic
    d = length(edge_nb1) - length(edge_nb2);
    if d > 0
        edge_nb2 = [edge_nb2; zeros([abs(d),1])];
    elseif d < 0
        edge_nb1 = [edge_nb1; zeros([abs(d),1])];
    end

    %%% Find edges with both nodes satisfying conditions
    % TODO: edge3plus is just a boolean with indices referencing
    % nb( prox_node_idx( prox_edges(:,1) ) )
    % This needs to be connected back to the original indices of nb. Then,
    % we need to extract nb( indices(edge3plus) ).
    edge3plus = (edge_nb1 > 2) & (edge_nb2 > 2);

    %%% Find edges (arrays of node indices) corresponding to edge3plus
    % edge3plus ([m, 1] logical) -- True/False for each entry in ref1, ref2
    % ref1, ref2 ([m, 1] double) -- proximal node indices
    % prox_edges ([n, 2] double) -- edges within proximity of bif3 node
    %       with values corresponding to indices in prox_node_idx
    
    %% Identify unique loops if there are any
    % TODO:
    % - how does this identify unique loops?
    % - why is nLoops calculated in this manner?
    nLoops = size(prox_edges,1) - length(prox_node_idx) + 1;
    if nLoops>0
        %% Find independent cycles of connected simple graph
        % Should prox_edges be replaced by the vertices of edge3plus?
        

        % Find basis of cycles of graph
        c = grCycleBasis( prox_edges );
        
        %% find loops with a unique 3+ edge
        % TODO: incompatible sizes. Should edge3plus be replaced with
        % proximal edges?
        lst = find( (sum(c,2)'.*edge3plus')==1 );
        
        % remove no more than one of these edges
        % from each loop and no nodes will be abandoned
        %
        % COULD THIS LEAVE A NODE WITH ONE EDGE IF I INADVERTANTLY REMOVE
        % EDGES FRO A GIVEN NODE TO REDUCE nB TO 1. I SHOULD CHECK THAT
        % nB>2 BEFORE REMOVING AN EDGE... EASY CHECK
        %
        if ~isempty(lst)
            % What is the purpose of lst3?
            lst3 = [];
            for iLoop = 1:nLoops
                % What is the purpose of lst2?
                lst2 = find(c(lst,iLoop).*edge3plus(lst)==1);
                if ~isempty(lst2)
                    % What is the purpose of iE?
                    iE = lst(lst2(1));
                    if nb(prox_node_idx(prox_edges(iE,1))) > 2 && nb(prox_node_idx(prox_edges(iE,2))) > 2
                        nb(prox_node_idx(prox_edges(iE,1)))...
                            = nb(prox_node_idx(prox_edges(iE,1))) - 1;
                        nb(prox_node_idx(prox_edges(iE,2)))...
                            = nb(prox_node_idx(prox_edges(iE,2))) - 1;
                        lst3(end+1) = iE;
                    end
                end
            end
            edges(prox_edge_idx(lst3),:) = [];
        else
            n1 = n1 + 1;
            n2 = n2 + nLoops;
        end
                
    end
end

% Update the Graph struct
Data.Graph.nB = nB;
Data.Graph.edges = edges;

% User Interface
close(hwait)
fprintf('No loops removed for %d nodes for a total of %d loops skipped',n1,n2);

% set(handles.textNumEdges,'string',sprintf('%d edges',size(Data.Graph.edges,1)))
% set(handles.uipanelGraph,'title',sprintf('Graph (%d nodes)',size(Data.Graph.nodes,1)));
% 
% set(handles.pushbuttonImageGraph,'enable','on')