%% Test regraphNodes with list of nodes
%{
Purpose:
- Etienne has data stored as .NII files
- this script will convert his .NII to .mat
%}
clear; close all; clc;

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
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));

%% Initialize data paths for dataset with loops