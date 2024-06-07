%% Test script for combining segmentations and graphing the results.
% Author: Mack Hyman
% Date Created: Dec. 13, 2023
%
%{
This script will iterate over a range of sigma values, which correspond to
the sigma value of the frangi filter. It will combine the segmentations
from each of these segmentations into a single matrix. Then, it will apply
a minimum probability threshold and then save it to the specified
directory.
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
subid = {'AD_20832', 'AD_20969','AD_21354', 'AD_21424',...
         'CTE_6489','CTE_6912','CTE_7019','CTE_7126','CTE_8572',...
         'NC_6839','NC_6974','NC_8653','NC_21499','NC_8095'};
subid = {'NC_8095'};

% Volume filename
volname = 'ref_4ds_norm_inv.tif';

% Sigma subdirectories ctontaining segmentation TIF files
sigma = {'gsigma_1-3-5_gsize_5-13-21','gsigma_2-3-4_gsize_9-13-17',...
        'gsigma_3-5-7_gsize_13-21-29','gsigma_5-7-9_gsize_21-29-37',...
        'gsigma_7--9-11_gsize_29-37-45'};

% Combined segmentation output folder name
dirout = 'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11';

%% Initialize struct for storing filepaths
ov = struct();
for ii = 1:length(subid)
    % Subject ID string
    ov(ii).subid = subid{ii};
    % Sigma arrays from Frangi
    ov(ii).s(1).sigma = sigma{1};
    ov(ii).s(2).sigma = sigma{2};
    ov(ii).s(3).sigma = sigma{3};
    ov(ii).s(4).sigma = sigma{4};
    ov(ii).s(5).sigma = sigma{5};
end

%%% Call function to create full file paths
ov = make_fpaths(dpath, subdir, ov, volname, dirout);

%% Threshold the probability maps and overlay segmentations
% Minimum threshold
th = 0.18;
% Subdirectory
subdir = '/p18/';
% Threshold probability maps and combine
combine_pmats(ov,th,subdir);

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
            ov(ii).f(j).fpath = fullfile(basepath, ...
                ov(ii).s(j).sigma, 'probability_map.mat');
        end
    end
end