%% Main file for calling converting segmentation to graph
% Author: Mack Hyman
% Date Created: March 12, 2024
%
% Detailed Description
%{
This script performs the following:
- convert segmentation to skeleton
- convert skeleton to graph
- remove loops from graph
%}
clear; clc; close all;

%% Convert 