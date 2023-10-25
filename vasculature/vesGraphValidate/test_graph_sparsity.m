%% Test script for the function graph_sparsity
% PURPOSE: Use synthetic data to validate the ability of the graph_sparsity
% to identify sparse cycles.
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

%% Call graph_sparsity
sp = graph_sparsity(edges);

%%% Plot graph
% Convert to graph and plot
g = graph(edges(:,1), edges(:,2));
visualize_graph(nodes, edges, 'Synthetic Test Data',[]);


%% Compare output to ground truth
% The synthetic data contain five cycles. The first four and last are
% sparse, and the last is not sparse. The output of "graph_sparsity" should
% match the array below.
ground_truth = [1; 1; 1; 1; 0; 1];
assert(all(ground_truth == sp), ['The graph_sparsity function did not return the' ...
    'correct answer'])
