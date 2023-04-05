%% Script for debugging the zero indexing issue in nodeGrps_vesSegment
% This script was created for issue #2 on the GitHub page:
% https://github.com/BUNPC/psoct_vessel_graphing/issues/2
clear; clc; close all;

%% Load data (these were generated by following normal procedure for GUI)
% These are the arrays nodeEdges (a,2) and nodePos (b,3)
load('nodeGrps_vesSegment_test_data.mat');

%% Call nodeGrps_vesSegment
nodeGrps_vesSegment(nodePos, nodeEdges);
