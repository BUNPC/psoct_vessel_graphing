%% Test script for the function rm_loop_edge
% PURPOSE: Use synthetic data to validate the ability of rm_loop_edge
% to remove one edge from a sparse loop.
% Author: Mack Hyman (Oct. 2023)

clear; clc; close all;
%% Generate synthetic data 
%%% SPARSE: Square w/ 4 nodes + 4 connecting segments at each node
s = [1 2 2 3 3 5 5 7];
t = [2 3 7 4 5 6 7 8];
e = vertcat(s,t);
e = e';
% Define node positions
n = [0 0; 2 0; 4 2; 4 4; 6 0; 8 0; 4 -2; 4 -4;];
% Add zeros for z position
n = horzcat(n, zeros(size(n,1),1));
% Create graph and visualize for verification
g = graph(s,t);
% visualize_graph(n, e, 'Sparse Square',[]);

%%% SPARSE: Square w/ diagonal + 4 connecting segments
s2 = s + 8;
t2 = t + 8;
e2 = vertcat(s2,t2);
e2 = e2';
e2 = [e2; 11, 15];
% Define node positions
n2 = n + [0,-10,0];

%%% Concatenate nodes and edges
etot = vertcat(e, e2);
ntot = vertcat(n, n2);

%%% NOT SPARSE: Square w/ extra + 4 connecting segments at each node
e3 = [17 18; 18 26; 26 19; 19 20; 19 27; 27 21; 21 22; 21 28; 28 23;...
    23 24; 23 25; 25 18];
% Define node positions
n3 = n(1:8,:) + [12,0,0];
n3 = vertcat(n3, n3(7,:) + [-1, 1, 0], n3(2,:) + [1, 1, 0],...
    n3(2,:) + [3, 1, 0], n3(7,:) + [1, 1, 0]);

%%% STRAIGHT LINE:
% s1 = [1 2 3 4 5 6];
% t1 = [2 3 4 5 6 7];
% % positions
% n1 = [1 1; 2 2; 3 3; 4 4; 5 5; 6 6; 7 7];

%%% Plot Final Results
% Concatenate all nodes/edges
etot = vertcat(etot, e3);
ntot = vertcat(ntot, n3);
% Convert to graph and plot
g = graph(etot(:,1), etot(:,2));
visualize_graph(ntot, etot, 'Synthetic Data',[]);

%% Identify sparse loops
sp = graph_sparsity(etot);

%% Remove longest edge 
% Find cycles in graph
c = allcycles(g);
% Convert sparsity array to boolean
sp = boolean(sp);
% Extract just the node indices from sparse cycles
csp = c;
csp(~sp) = [];


% Remove the longest edge from each
e_rm = rm_loop_edge(ntot, etot, sp, csp);

% Plot results
visualize_graph(ntot, e_rm, {'Synthetic Data','After Longest Edge Removal'},[]);