%% Test script for combining segmentations and graphing the results.
% Author: Mack Hyman
% Date Created: Dec. 13, 2023
%
%{
Detailed Description:
This script will pass several segmentations to the function
"combine_segs_then_graph" which will then combine the segmentations,
initialize the graph, remove loops, and save the output.

To Do:
- Finish function "make_fpaths"
- finish function "combine_segs_then_graph"
%}
clear; clc; close all;

%% Add top-level directory of code repository to path
% This allows Matlab to find the functions in the project folders

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

%% Initialization parameters for all subjects

% Struct for initializing filepaths
ov = struct();

if ispc
    % Top-level path to data
    dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\';
    % Subdirectory to data
    subdir = '\dist_corrected\volume\';

    %%% Subject IDs
    ov(1).subid = 'AD_10382';
    ov(2).subid = 'AD_20832';
    
    %%% Sigma subfolders for each Frangi sigma array
    ov(1).s(1).sigma = 'gsigma_1-3-5_gsize_5-13-21';
    ov(1).s(2).sigma = 'gsigma_5-7-9_gsize_21-29-37';
    ov(1).s(3).sigma = 'gsigma_7--9-11_gsize_29-37-45';
    % For now, copy these into the second subject
    ov(2).s = ov(1).s;

    %%% TIF filename (w/ probability map threshold) for each subject
    ov(1).f(1).fname = 'ref_4ds_norm_inv_segment_pmin_0.23.tif';
    ov(1).f(2).fname = 'ref_4ds_norm_inv_segment_pmin_0.23.tif';
    ov(1).f(3).fname = 'ref_4ds_norm_inv_segment_pmin_0.23.tif';
    ov(2).f(1).fname = 'ref_4ds_norm_inv_segment_pmin_0.23.tif';
    ov(2).f(2).fname = 'ref_4ds_norm_inv_segment_pmin_0.23.tif';
    ov(2).f(3).fname = 'ref_4ds_norm_inv_segment_pmin_0.23.tif';
    
    %%% Call function to create full file path
    ov = make_fpaths(dpath, subdir, ov);
else
    error('Just a local test script for debugging. Jog on, mate.')
end

%% Call function to overlay

%%% String for main directory containing all data









%% Function to create filepaths
function [overlay] = make_fpaths(dpath, subdir, ov)
    %%% Iterate over each subject
    for ii = 1:size(ov, 2)
        % Create base path to data
        basepath = fullfile(dpath, subdir);
      
        %%% Iterate over sigma arrays
        for j = 1:length(ov(ii).s)
            % Create filepath for each sigma
            % TODO: update output struct (ov_final(?,?))
            overlay(ii).f(j).fpath = fullfile(basepath, ov(ii).s(j).sigma, ov(ii).f(j).fname);
        end
    end
end




























%% Have a wonderful day, and go f