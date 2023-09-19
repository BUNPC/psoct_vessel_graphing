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
% - find nodes only belonging to loops
%       - will need to update "remove_loops" after this update

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
%%% Data with fewer nested loops (probability threshold = 0.23)
seg_name = 'ref_4ds_norm_inv_crop2_segment_pmin_0.23.tif';
% Graphed data after Gaussian filtering segmentation
gdata = 'ref_4ds_norm_inv_crop2_segment_pmin_0.23_mask40_graph_data.mat';
%%% Data with nested loops (probability threshold = 0.21)
% seg_name = 'ref_4ds_norm_inv_crop2_segment_pmin_0.21.tif';
% gdata = 'ref_4ds_norm_inv_crop2_segment_pmin_0.21_mask_40_graph_data.mat';

%% Load PSOCT graph
Data = load(fullfile(dpath, subid, subdir, sigdir, gdata), 'Data');
Data = Data.Data;

%%% Create new variable to match format of function
% im.angio = tmp.angio;
im.nodes = Data.Graph.nodes;
im.edges = Data.Graph.edges;
im.segn = Data.Graph.segInfo.nodeSegN;

%%% Load volumetric information and set threshold
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

%%% Load segmentation stack
% Import volume
seg = TIFF2MAT(fullfile(dpath, subid, subdir, sigdir, seg_name));
volshow(seg);

%%% Parameters to convert graph to 3D skeleton
%  sz =  the size of output volume. The order is [ y x z ]
sz = size(seg);
%  res = resolution of the nifti image [y,x,z] centimeters
res = [0.0012, 0.0012, 0.0015];
%  ds_flag = 1 or 0. (1=downsampled. 0=no-downsampled)
ds_flag = 1;
%  save_flag = to save the skeleton or not
save_flag = 0;


%% Call function to remove loops (regraph + move to mean)

%%% [x,y,z] limits for visualizing graph
% Cropped limits of branch point
% xlims = 100:180; ylims = 130:220; zlims = 60:80;
% Cropped limits of graph with loops
% xlims = 1:250; ylims = 1:30; zlims = 120:140;
% Cropped limits of tiered loops
xlims = 150:300; ylims = 120:300; zlims = 10:50;

%%% Visualize graph prior to processing
graph_vis(im.nodes, im.edges, 'Graph Before Processing');
% xlim(xlims); ylim(ylims); zlim(zlims);

%%% Call remove loops function
% Minimum voxel intensity threshold for moving to mean
vmin = 0.5;
[n, e] = remove_loops(im.nodes, im.edges, vol, vmin);

%%% Plot results
t_str = strcat("Delta = 2. ",'Voxel Intensity Threshold = ',num2str(vmin));
g_title = {'Regraphed & Moved to Mean', t_str};
graph_vis(n, e, g_title)
% xlim(xlims); ylim(ylims); zlim(zlims);

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

%% Overlay segmentation w/ skeleton (from centered graph)
%{
%%% Convert graph to 3D skeleton
centered_fout = strcat(seg_name, 'centered.nii');
[centered_skel] = sk3D(sz, im_centered, centered_fout,...
                        res, ds_flag, save_flag);

%%% Overlay in figure
view_panel = uifigure(figure,'Name',"Centered");
close
v = viewer3d(view_panel);
v.BackgroundColor = 'w';
v.BackgroundGradient = 'off';
h = volshow(centered_skel,'Parent',v);
h.Parent.BackgroundColor = 'w';
h.OverlayData = seg;
h.OverlayAlphamap = 0.3;

%%% Create a .gif of rotating the overlaid volumes
%{
gif_filename = 'regraph_mv2mean_regraph_center.gif';
gif_fullfile = fullfile(dpath, subid, subdir, sigdir,...
    'ref_4ds_norm_inv_crop2_segment_pmin_0.23_gifs', gif_filename);
volgif(h, centered_skel, gif_fullfile);
%}
%}

%% Convert uncentered graph to 3D skeleton
% Output filename
uncentered_fout = strcat(seg_name, 'uncentered.nii');
% Convert
[uncentered_skel] = sk3D(sz, im_re, uncentered_fout,...
                        res, ds_flag, save_flag);

%% Overlay uncentered & centered
%{
view_panel = uifigure(figure,'Name',"Uncentered"); close;
v = viewer3d(view_panel);
v.BackgroundColor = 'w';
v.BackgroundGradient = 'off';
h = volshow(uncentered_skel,'Parent',v);
h.Parent.BackgroundColor = 'w';
h.OverlayData = centered_skel;
h.OverlayAlphamap = 0.3;
%}

%% Overlay segmentation w/ processed (but uncentered) skeleton
% GOAL: ensure we can remove loops without modifying branch points. This
% was verified visually that it does not remove branch points.

% Cropped limits of branch point
% xlims = 100:180; ylims = 130:220; zlims = 60:80;
% Cropped limits of graph with loops
xlims = 1:250; ylims = 1:30; zlims = 120:140;
% Crop the segmentation into a subvolume 
seg_crop = seg(xlims, ylims, zlims);
proc_skel_crop = uncentered_skel(xlims, ylims, zlims);

view_panel = uifigure(figure,'Name',"Segmentation + Uncentered"); close;
v = viewer3d(view_panel);
v.BackgroundColor = 'w';
v.BackgroundGradient = 'off';
h = volshow(proc_skel_crop,'Parent',v);
h.Parent.BackgroundColor = 'w';
h.OverlayData = seg_crop;
h.OverlayAlphamap = 0.1;

%% Overlay segmentation + raw graph
%%% Convert graph to 3D skeleton
% Output filename
nonproc_fout = strcat(seg_name, 'raw_graph_skel.nii');
% Convert unprocessed skeleton
[raw_skel] = sk3D(sz, im, nonproc_fout,...
                        res, ds_flag, save_flag);
% Crop the raw skeleton
raw_skel_crop = raw_skel(xlims, ylims, zlims);

%%% Overlay segmentation + raw graph
view_panel = uifigure(figure,'Name',"Unprocessed"); close;
v = viewer3d(view_panel);
v.BackgroundColor = 'w';
v.BackgroundGradient = 'off';
h = volshow(raw_skel_crop,'Parent',v);
h.Parent.BackgroundColor = 'w';
h.OverlayData = seg_crop;
h.OverlayAlphamap = 0.1;

%% Overlay skeletons and segmentation from cropped volume
% Cropped limits of tiered loops
xlims = 150:300; ylims = 120:300; zlims = 10:50;
% Crop the segmentation into a subvolume 
seg_crop = seg(xlims, ylims, zlims);

%%% Crop unprocessed
raw_skel_crop = raw_skel(xlims, ylims, zlims);
view_panel = uifigure(figure,'Name',"Unprocessed + Segmentation"); close;
v = viewer3d(view_panel);
v.BackgroundColor = 'w';
v.BackgroundGradient = 'off';
h = volshow(raw_skel_crop,'Parent',v);
h.Parent.BackgroundColor = 'w';
h.OverlayData = seg_crop;
h.OverlayAlphamap = 0.1;

%%% Crop processed (without centering)
proc_skel_crop = uncentered_skel(xlims, ylims, zlims);
volshow(proc_skel_crop);
view_panel = uifigure(figure,'Name',"Processed + Segmentation"); close;
v = viewer3d(view_panel);
v.BackgroundColor = 'w';
v.BackgroundGradient = 'off';
h = volshow(proc_skel_crop,'Parent',v);
h.Parent.BackgroundColor = 'w';
h.OverlayData = seg_crop;
h.OverlayAlphamap = 0.1;

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

%% Create .gif of volume rotating
function volgif(h, volume, filename)
% INPUTS:
%   h (handle): handle to the volshow output
%   volume (matrix): matrix of the volume in volshow
%   filename (string): full output file path + filename
% OUTPUTS:
%   the .gif is saved to the filename specified in "fname"

viewer = h.Parent;
hFig = viewer.Parent;
drawnow

% Aim camera at center of volume
sz = size(volume);
center = sz/2 + 0.5;
viewer.CameraTarget = center;

% Specify number of frames in animation
nframes = 20;
vec = linspace(0,2*pi,nframes)';
dist = sqrt(sz(1)^2 + sz(2)^2 + sz(3)^2);
myPosition = center + ([cos(vec) sin(vec) ones(size(vec))]*dist);

% Rotate camera, update display, write to gif
for idx = 1:length(vec)
    % Update the current view
    viewer.CameraPosition = myPosition(idx,:);
    % Capture the image using the getframe function
    I = getframe(hFig);
    [indI,cm] = rgb2ind(I.cdata,256);
    % Write the frame to the GIF file
    if idx==1
        % Do nothing. The first frame displays only the viewer, not the
        % volume.
    elseif idx == 2
        imwrite(indI,cm,filename,"gif",Loopcount=inf,DelayTime=0)
    else
        imwrite(indI,cm,filename,"gif",WriteMode="append",DelayTime=0)
    end
end

end

%% Regraph and move to mean until all loops removed

function [nodes, edges] = remove_loops(nodes, edges, angio, v_min)
%rm_loops Remove loops in graph.
%   Outline:
%       - Use graph function "allcycles" to find loops
%       - While loops exist:
%           - convert from graph back to nodes/edges
%           - regraph (downsample nodes within search radius)
%           - move to mean (collapse loops)
%
%   INPUTS:
%       nodes ([n,3] array): node locations
%       edges ([m,2] array): edges connecting each node
%       im.angio (double matrix): PS-OCT intensity volume (vessels are
%               bright)
%       v_min (double): minimum voxel intensity threshold. The new voxel
%               position will only be reassigned if the voxel intensity of
%               the new node position is >= v_min.
%   OUTPUTS:
%       n ([n,3] array): node locations
%       e ([m,2] array): edges connecting each node

%% Convert from nodes + edges into Matlab graph
% Create standard Matlab graph
g = graph(edges(:,1), edges(:,2));

% Detect loops in graph
cycles = allcycles(g);

%% While loops exist: regraph + move to mean
while ~isempty(cycles)
    %%% Regraph (downsample) to remove collapsed loops
    % Initialized validated nodes vector so that regraph code runs
    validated_nodes = zeros(size(nodes,1),1);
    % Search delta for the x,y,z directions (units = voxels)
    delta = 2;
    % Call function to regraph
    [nodes, edges, ~, ~] =...
        regraphNodes_new(nodes, edges, validated_nodes, delta);
    
    %%% Move to mean (collapse loops)
    % Create copy of regraphed graph
    im_mv = struct();
    im_mv.nodes = nodes;
    im_mv.edges = edges;
    im_mv.angio = angio;
    % Perform move to mean five times to collapse nodes. This number was
    % determined emprically from testing a data subset.
    for j=1:5
        im_mv = mv_to_mean(im_mv, v_min);
    end
    
    %%% Regraph to smooth collapse loops
    [nodes, edges, ~, ~] =...
        regraphNodes_new(im_mv.nodes, im_mv.edges, validated_nodes, delta);
    
    %%% Detect loops in graph
    g = graph(edges(:,1), edges(:,2));
    cycles = allcycles(g);
end

end

















