%% Test script - yaxis_scatter
% This requires the angio, node positions, edges, threshold for
% segmentation, and the voxel size.

clear; clc; close all;

%% Add top-level directory of code repository to path
% This allows Matlab to find the functions in the project folders
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

%% Define directory paths

% Path to top-level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
   
% Volume and graph directories
voldir = '/dist_corrected/volume/ref/';
maskdir = 'dist_corrected/volume/ref/masks/';
graphdir = ['/dist_corrected/volume/combined_segs/' ...
            'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];

% Metrics output path
mpath = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
    'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];


% IDs of each subject
subids = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_7126',...
         'NC_6839',  'NC_6974',  'NC_8653',  'NC_21499', 'NC_8095'};

%% Organize metrics into disease and control group
%%% Initialize metrics and regions
% Load the struct with average metrics for each region
load(fullfile(mpath,'metrics.mat'));
% Average Metric Parameters
params = {'length_density','branch_density','fraction_volume',...
    'tortuosity'};
% Region of brain (excluding ratios)
regions = {'tiss','gyri','sulci','gm','wm','gm_sulci','wm_sulci',...
            'gm_gyri','wm_gyri'};

%%% Iterate over each tissue region
for ii = 1:length(regions)
    %%% Organize the average metric for each subject
    % This includes the array of tortuosity values for each subject
    for j = 1:length(params)
        % Retrieve the average metric values for this region/parameter from
        % all AD, CTE, and NC subjects. Concatenate into arrays.
        [ad, cte, nc] = organize_metrics(metrics, subids,...
                                         regions{ii}, params{j});
        % Save parameter by the group to fascilitate statistical analyses
        % in a later step.
        metrics.(regions{ii}).(params{j}).ad = ad;
        metrics.(regions{ii}).(params{j}).cte = cte;
        metrics.(regions{ii}).(params{j}).nc = nc;
    end
end

%% Generate y-axis scatter plot for one of the metrics

% Colors of plots
c = [0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.8500 0.3250 0.0980];

% Iterate over vascular metrics
for ii = 1:length(params)
    % Create array for each region
    for j = 1:length(regions)
        % Retrieve array for each group
        ad = metrics.(regions{j}).(params{ii}).ad;
        cte = metrics.(regions{j}).(params{ii}).cte;
        nc = metrics.(regions{j}).(params{ii}).nc;
        % scattter plot
        yaxis_scatter(c)

    end
end













