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
seg_name = 'ref_4ds_norm_inv_crop2_segment_pmin_0.23.tif';
% Graphd data after Gaussian filtering segmentation
gdata = 'ref_4ds_norm_inv_crop2_segment_pmin_0.23_mask40_graph_data.mat';

%% Load PSOCT graph
Data = load(fullfile(dpath, subid, subdir, sigdir, gdata), 'Data');
Data = Data.Data;

% Create new variable to match format of function
% im.angio = tmp.angio;
im.nodes = Data.Graph.nodes;
im.edges = Data.Graph.edges;
im.segn = Data.Graph.segInfo.nodeSegN;

% Set x,y,z limits for graph
xlims = [0, 30];
ylims = [0, 250];
zlims = [120, 140];
%% Load volumetric information and set threshold
% Import volume
vol = TIFF2MAT(fullfile(dpath, subid, subdir, vdata));
im.angio = vol;
% volshow(vol);

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

%% Load segmentation stack
% Import volume
seg = TIFF2MAT(fullfile(dpath, subid, subdir, sigdir, seg_name));
volshow(seg);

%% Visualize graph prior to processing
graph_vis(im.nodes, im.edges, 'Graph Before Processing');
xlim(xlims); ylim(ylims); zlim(zlims);
set(gca, 'FontSize', 25);
%% Downsample (regraph)
% Initialized validated nodes vector so that regraph code runs
validated_nodes = zeros(size(im.nodes,1),1);
% Search radius (units = voxels)
delta = 2;

% Call regraph (downsample)
[im_ds.segn, im_ds.nodes, im_ds.edges, ~, ~] =...
    regraphNodes_new(im.segn, im.nodes, im.edges, validated_nodes, delta);

% Add angio to new struct
im_ds.angio = im.angio;

% Visualize downsampled
gt1 = {'Regraphed', strcat("Delta = ", num2str(delta))};
graph_vis(im_ds.nodes, im_ds.edges, gt1);
xlim(xlims); ylim(ylims); zlim(zlims);
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
xlim(xlims); ylim(ylims); zlim(zlims);
%% Downsampling (regraph) to collapse loop
validated_nodes = zeros(size(im.nodes,1),1);
% Regraph
[im_re.segn, im_re.nodes, im_re.edges, ~, ~] =...
    regraphNodes_new(im_mv.segn, im_mv.nodes, im_mv.edges, validated_nodes, delta);
% Assign angio
im_re.angio = vol;

% Graph title
t_str = strcat("Regraph Delta = ", num2str(delta),...
    '. Voxel Intensity Threshold = ',num2str(th_norm(ii)));
% Visualize title
graph_vis(im_re.nodes, im_re.edges, {'Regraphed, Moved to Mean, Regraphed', t_str})
% xlim(xlims); ylim(ylims); zlim(zlims);
%% Centering (determine which centering function to use)
% There are several different centering functions. David recommended trying
% or reviewing each to determine which is most suitable.

%%% Replace the angio (PS-OCT volume) with the segmentation
% This is a debugging step. The PSOCT has poor contrast, so centering does
% not work. David recommended centering with the segmentation.
im_re.angio = seg;
% im_re.angio = vol;

% Normalize posM
centerStep1vox = 0;

% Visualization
visualize_flag = 0;

% Run centering function
im_centered = center_nodes_xyz(im_re, centerStep1vox, visualize_flag, im_thresh);

%% Graph after centering
% Visualize centered graph
graph_vis(im_centered.nodes, im_centered.edges, 'Graph After Centering')
xlim(xlims); ylim(ylims); zlim(zlims);
% xlim([200, 300]); ylim([150,300]); zlim([50, 150]);

%% Overlay centering with segmentation
% GOAL: ensure we can remove loops without modifying branch points
% - Verify it does not grossly modify/eliminate branch points
% - Determine if centering moves nodes towards center line.

%%% Convert results from graph to skeleton (save nifti)
% Output filename
centered_fout = strcat(segname, 'centered');
%  sz =  the size of output volume. The order is [ y x z ]

%  Graph  =  Graph struct that has node coordinate vector and edge vector 

%  res    =  resolution of the nifti image. 

%            i.e. [0.01, 0.01, 0.01] 0.01cm (10um) isotropics

%  downsampled_flag = 1=downsampled; 0=no-downsampled

%  save_flag = to save the skeleton or not

[~] = sk3D(sz,Graph,centered_output,res,downsampled_flag,save_flag);

%% Overlay centering results with uncentered

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
p.LineWidth = 4;

% initialize camera view
view(3);

% Labels, title, fontsize, grid
title(title_str); xlabel('x'); ylabel('y'); zlabel('z')
set(gca, 'FontSize', 25);
grid on;
end

%% Separate 