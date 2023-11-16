%% Generate synthetic graph data for debugging
% PURPOSE: Create synthetic data to validate graph_sparsity, rm_loop_edge
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
% Concatenate nodes and edges
etot = vertcat(e, e2);
ntot = vertcat(n, n2);

%%% NOT SPARSE: Square w/ extra + 4 connecting segments at each node
e3 = [17 18; 18 26; 26 19; 19 20; 19 27; 27 21; 21 22; 21 28; 28 23;...
    23 24; 23 25; 25 18];
% Define node positions
n3 = n(1:8,:) + [12,0,0];
n3 = vertcat(n3, n3(7,:) + [-1, 1, 0], n3(2,:) + [1, 1, 0],...
    n3(2,:) + [3, 1, 0], n3(7,:) + [1, 1, 0]);
% Concatenate all nodes/edges
etot = vertcat(etot, e3);
ntot = vertcat(ntot, n3);

%%% NOT SPARSE: assymetrical rectangle with one longest edge
% Edges
e = [29 30; 30 31; 31 32; 31 38; 38 39; 39 40; 38 35; 32 35; 32 33;...
    33 34; 35 36; 36 37];
% Node positions
n = [13 10; 14 10; 15 10; 18 8; 17 6; 17 5; 19 10; 20 10; 21 10;...
    17 12; 17 13; 17 14];
n(:,2) = n(:,2) .* -1;
n = horzcat(n, zeros(length(n),1));
% Concatenate to overall list
etot = vertcat(etot, e);
ntot = vertcat(ntot, n);

%%% STRAIGHT LINE:
% s1 = [1 2 3 4 5 6];
% t1 = [2 3 4 5 6 7];
% % positions
% n1 = [1 1; 2 2; 3 3; 4 4; 5 5; 6 6; 7 7];

%%% Plot Final Results
% Convert to graph and plot
g = graph(etot(:,1), etot(:,2));
visualize_graph(ntot, etot, 'Synthetic Test Data',[]);


%% Save output as .MAT file
% Rename variables for usability
edges = etot;
nodes = ntot;
save('synthetic_graph.mat', 'edges', 'nodes');
