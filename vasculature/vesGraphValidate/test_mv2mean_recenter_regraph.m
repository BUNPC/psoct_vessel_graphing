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

%%% Programatically compared uncentered w/ centered
nodes_cent = im_centered.nodes;
nodes_uncent = im_re.nodes;
cent_comp = nodes_cent == nodes_uncent;

%% Graph after centering
% Visualize centered graph
graph_vis(im_centered.nodes, im_centered.edges, 'Graph After Centering')
% xlim(xlims); ylim(ylims); zlim(zlims);
% xlim([200, 300]); ylim([150,300]); zlim([50, 150]);

%% Overlay segmentation w/ skeleton (from centered graph)
% GOAL: ensure we can remove loops without modifying branch points
% - Verify it does not grossly modify/eliminate branch points
% - Determine if centering moves nodes towards center line.

%%% Parameters to convert graph to 3D skeleton
% Output filename
centered_fout = strcat(seg_name, 'centered.nii');
%  sz =  the size of output volume. The order is [ y x z ]
sz = size(seg);
%  res = resolution of the nifti image [y,x,z] centimeters
res = [0.0012, 0.0012, 0.0015];
%  ds_flag = 1 or 0. (1=downsampled. 0=no-downsampled)
ds_flag = 1;
%  save_flag = to save the skeleton or not
save_flag = 0;

%%% Convert graph to 3D skeleton
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

%% Overlay uncentered & centered
%%% Convert graph to 3D skeleton
% Output filename
uncentered_fout = strcat(seg_name, 'uncentered.nii');
% Convert
[uncentered_skel] = sk3D(sz, im_centered, uncentered_fout,...
                        res, ds_flag, save_flag);

%%% Overlay in figure
h = volshow(uncentered_skel);
h.OverlayData = centered_skel;
h.OverlayAlphamap = 0.3;

%% Overlay segmentation w/ skeleton (from uncentered graph)
h = volshow(uncentered_skel);
h.OverlayData = seg;
h.OverlayAlphamap = 0.3;

%% Overlay segmentation + unprocessed graph
%%% Convert graph to 3D skeleton
% Output filename
nonproc_fout = strcat(seg_name, 'raw_graph_skel.nii');
% Convert
[raw_skel] = sk3D(sz, im_centered, nonproc_fout,...
                        res, ds_flag, save_flag);

%%% Overlay segmentation + raw graph
view_panel = uifigure(figure,'Name',"Unprocessed");
close
v = viewer3d(view_panel);
v.BackgroundColor = 'w';
v.BackgroundGradient = 'off';
h = volshow(raw_skel,'Parent',v);
h.Parent.BackgroundColor = 'w';
h.OverlayData = seg;
h.OverlayAlphamap = 0.3;

%% Overlay cropped raw skeleton
raw_skel_crop = raw_skel(100:300, 100:300,100:169);
volshow(raw_skel_crop);

proc_skel_crop = uncentered_skel(100:300, 100:300,100:169);
volshow(proc_skel_crop);

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



