%% Test script for the CenterNodesXYZ function from David Boas
% David wrote code for moving nodes towards the center line of a segment.
% The original function was integrated into a GUI, but there is no
% documentation for the GUI. Therefore, this script will test a standalone
% version of the function.
%
% This script is currently removing loops via the following steps:
% Regraph (delta = 2)
% Move to mean (normalized voxel intensity threshold = 0.5)
% Regraph (delta = 2)
% 
% The main issue with this script is that the downsampling (regraph) is
% combining disparate segments. This requires updating the regraph to only
% downsample nodes along the same edge.

% TODO:
% - modify regraph to only downsample nodes in same segment

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

%% Setup filename
% Paths to files
dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\';
subid = 'NC_6839';
subdir = '\dist_corrected\volume\';
sigdir = 'gsigma_1-3-5_gsize_5-13-21\';
vdata = 'ref_4ds_norm_inv_crop2.tif';
% Graphd data after Gaussian filtering segmentation
gdata = 'ref_4ds_norm_inv_crop2_segment_pmin_0.23_mask40_graph_data.mat';

%% Load PSOCT graph
tmp = load(fullfile(dpath, subid, subdir, sigdir, gdata), 'Data');
tmp = tmp.Data;

% Create new variable to match format of function
% im.angio = tmp.angio;
im.nodes = tmp.Graph.nodes;
im.edges = tmp.Graph.edges;

%% Load volumetric information and set threshold

% Import volume
vol = TIFF2MAT(fullfile(dpath, subid, subdir, vdata));
im.angio = vol;

% Imaging voxel intensity threshold (normalized [0,1])
im_thresh = 0.75;

% Check for data type
if strcmp(class(vol),'uint16')
    r = 65535;
    im_thresh = im_thresh .* r;
elseif strcmp(class(vol),'uint8')
    r = 255;
    im_thresh = im_thresh .* r;
end

%% Visualize graph prior to processing
graph_vis(im.nodes, im.edges, 'Graph Before Processing')

%% Downsample (regraph)
% Initialized validated nodes vector so that regraph code runs
validated_nodes = zeros(size(im.nodes,1),1);
% Search radius (units = voxels)
delta = 2;

% Call regraph (downsample)
[im_ds.nodes, im_ds.edges, ~, ~] =...
    regraphNodes_new(im.nodes, im.edges, validated_nodes, delta);

% Add angio to new struct
im_ds.angio = im.angio;

% Visualize downsampled
gt1 = {'Regraphed', strcat("Delta = ", num2str(delta))};
graph_vis(im_ds.nodes, im_ds.edges, gt1);

%% Test move to mean and thresholding
th_norm = [0.25, 0.5, 0.75, 0.9, 0.95];
tarray = th_norm .* r;

for ii = 2:2
    im_mv = im_ds;
    for j=1:5
        im_mv = mv_to_mean(im_mv, tarray(ii));
    end
    %%% Visualize centered graph
    % Title for graph
    t_str = strcat("Delta = ", num2str(delta),...
        '. Voxel Intensity Threshold = ',num2str(th_norm(ii)));
    g_title = {'Regraphed & Moved to Mean', t_str};
    graph_vis(im_mv.nodes, im_mv.edges, g_title)
end

%% Downsampling (regraph) to collapse loop
validated_nodes = zeros(size(im.nodes,1),1);
[im_re.nodes, im_re.edges, ~, ~] =...
    regraphNodes_new(im_mv.nodes, im_mv.edges, validated_nodes, delta);

t_str = strcat("Delta = ", num2str(delta),...
    '. Voxel Intensity Threshold = ',num2str(th_norm(ii)),...
    ". Delta = ", num2str(delta));

graph_vis(im_re.nodes, im_re.edges, {'Regraphed, Moved to Mean, Regraphed', t_str})

%% Centering (determine which centering function to use)
% There are several different centering functions. David recommended trying
% or reviewing each to determine which is most suitable.

% Unsure what this variable does
centerStep1vox = 0;
% Visualization
visualize_flag = 0;
% Run centering function
im_centered = center_nodes_xyz(im, centerStep1vox, visualize_flag, im_thresh);

%%% Compare results
% Visualize uncentered graph
graph_vis(im.nodes, im.edges, 'Graph Before Centering')
% Visualize centered graph
graph_vis(im_centered.nodes, im_centered.edges, 'Graph After Centering')

%% Visualize graph
function graph_vis(nodes, edges, title_str)
% Copy edges into standard format

s = edges(:,1); % source node
t = edges(:,2); % target node

% Create standard Matlab graph
g = graph(s, t);

% Plot graph before removing loops
figure;
p = plot(g, 'XData', nodes(:,1), 'YData', nodes(:,2), 'ZData', nodes(:,3));
% Set nodes red
p.NodeColor = 'red';
% Set line width of edges
p.LineWidth = 2;
% Labels and title
title(title_str); xlabel('x'); ylabel('y'); zlabel('z')
% initialize camera view
view(3);
end