%% Test script for the function rm_loop_edge
% PURPOSE: Use synthetic data to validate the ability of rm_loop_edge
% to remove one edge from a sparse loop.
% Author: Mack Hyman (Oct. 2023)

clear; clc; close all;

%% Add top-level directory of code repository to path
% Start in current directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Truncate path to reach top-level directory (psoct_vessel_graphing)
topdir = mydir(1:idcs(end));
addpath(genpath(topdir));

%% Load synthetic graph data (edges, nodes)
fname = fullfile(pwd, 'synthetic_graph.mat');
load(fname);
visualize_graph(nodes, edges, {'Synthetic Data','Before Longest Edge Removal'},[]);

%% Identify sparse loops
sp = graph_sparsity(edges);
% Convert sparsity array to boolean
sp = boolean(sp);

%% Remove longest edge 
% Generate graph
g = graph(edges(:,1), edges(:,2));
% Find cycles in graph
[cnodes, ~] = allcycles(g);
% Keep node indices from sparse cycles
cnodes(~sp) = [];

% Remove the longest edge from each
e_rm = rm_loop_edge(nodes, edges, sp, cnodes);

%% Verify that all sparse cycles had edge removed
% Plot results
visualize_graph(nodes, e_rm, {'Synthetic Data','After Longest Edge Removed'},[]);

% Assert no sparse cycles
sp = graph_sparsity(e_rm);
assert(sp == 0, 'The rm_loop_edge function did not remove all sparse loops.')
