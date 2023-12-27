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
topdir = mydir(1:idcs(end));
addpath(genpath(topdir));

%% Initialize directory paths for all subjects

if ispc
    % Top-level path to data
    dpath = ['C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\' ...
        'test_data\Ann_Mckee_samples_10T\'];  
    % Subdirectory to volume 
    subdir = '\dist_corrected\volume\';
elseif isunix
    % Top-level path to data
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T';
    % Subdirectory to volume 
    subdir = '/dist_corrected/volume/';
end

% Subject IDs
subid = {'AD_10382', 'AD_20832', 'AD_20969',...
             'AD_21354', 'AD_21424',...
             'CTE_6489','CTE_6912',...
             'CTE_7019','CTE_8572','CTE_7126',...
             'NC_6047', 'NC_6839',...
             'NC_6974', 'NC_7597',...
             'NC_8095', 'NC_8653',...
             'NC_21499','NC_301181'};

% Volume filename
volname = 'ref_4ds_norm_inv.tif';

% Sigma subdirectories ctontaining segmentation TIF files
sigma = {'gsigma_3-5-7_gsize_13-21-29','gsigma_5-7-9_gsize_21-29-37',...
        'gsigma_7--9-11_gsize_29-37-45'};

% Combined segmentation output folder name
dirout = 'gsigma_3-5-7_5-7-9_7-9-11';

%% Initialize struct for storing filepaths
ov = struct();
for ii = 1:length(subid)
    % Subject ID string
    ov(ii).subid = subid{ii};
    % Sigma arrays from Frangi
    ov(ii).s(1).sigma = sigma{1};
    ov(ii).s(2).sigma = sigma{2};
    ov(ii).s(3).sigma = sigma{3};
    % Filename of segmentation
    ov(ii).f(1).fname = 'ref_4ds_norm_inv_segment_pmin_0.23.tif';
    ov(ii).f(2).fname = 'ref_4ds_norm_inv_segment_pmin_0.23.tif';
    ov(ii).f(3).fname = 'ref_4ds_norm_inv_segment_pmin_0.23.tif';
end

%%% Call function to create full file paths
ov = make_fpaths(dpath, subdir, ov, volname, dirout);

%% Call function to overlay
combine_segs(ov);

%% Function to create filepaths
function [ov] = make_fpaths(dpath, subdir, ov, volname, dirout)
    %%% Iterate over each subject
    for ii = 1:size(ov, 2)     
        % Create base path to data
        basepath = fullfile(dpath, ov(ii).subid, subdir);    

        % Add path to volume
        ov(ii).vol = fullfile(basepath, volname);

        % New directory path for segmentation output
        ov(ii).dirout = fullfile(basepath,'combined_segs', dirout);

        %%% Iterate over sigma arrays
        for j = 1:length(ov(ii).s)
            % Create filepath for each sigma
            % TODO: update output struct (ov_final(?,?))
            ov(ii).f(j).fpath = fullfile(basepath, ...
                ov(ii).s(j).sigma, ov(ii).f(j).fname);
        end
    end
end


%% Have a wonderful day, and go f