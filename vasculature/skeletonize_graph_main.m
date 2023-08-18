%% Main file for skeletonizing/graphing the vasculature
%
%{
This script performs the following:
- Calculate volume of masked PSOCT volume
- skeletonize/graph the segmentation
- calculate vascular metrics from graph

TODO:
- crop the tissue volume and rerun segmentation
- find optimal sigma for gaussian3d smoothing. This is currently bluring
together branches, which is unacceptable.
- attempt edge-preserving filtering (eg bilateral filtering)
https://www.mathworks.com/help/images/linear-filtering.html
- IDEA: 
    - edge preserve filter on each slice
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

%% Initialize data path for linux or personal machine (debugging)
% This section is for setting the directory path to your datasets.
% The code is written to create a filepath with the following structure:
% [dpath]\[subid]\[subdir]\[fname].[ext]

% Your local directory to the segmentation data
dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\';
% Subject IDs
subid = 'NC_6839';
subdir = '\dist_corrected\volume\gsigma_1-3-5_gsize_5-13-21\';
% Filenames to parse (test data)
vol_name = 'ref_4ds_norm_inv_crop2.tif';
angio_name = 'ref_4ds_norm_inv_crop2_segment_pmin_0.23_mask40.mat';
% Voxel dimensions (verify this with Dylan or another RA at Martinos)
vox_dim = [12, 12, 15];
% Volume of a single voxel (microns cubed)
vox_vol = vox_dim(1) .* vox_dim(2) .* vox_dim(3);
% Output filepath + filenames
fout = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\NC_6839\dist_corrected\volume\gsigma_1-3-5_gsize_5-13-21\';
angio_pre_fname = 'angio_pre.tif';
angio_post_fname = 'angio_post.tif';

%% Load the masked PSOCT volume
volpath = fullfile(dpath, subid, subdir);
vol_fname = fullfile(dpath, subid, '\dist_corrected\volume\', vol_name);
angio_fname = strcat(volpath, angio_name);

%%% Import volume
% Import tissue and convert .tif to .MAT
tissue = TIFF2MAT(vol_fname);
% Invert so that non-tissue voxels are zeros
tissue_inv = imcomplement(tissue);
% Calculate total number of non-zero voxels in tissue sample
tissue_logical = logical(tissue_inv);
voxels = sum(tissue_logical(:));
% Convert voxels to metric volume (cubic microns)
vol = voxels .* vox_vol;

%%% Import angio (segmentation)
angio = load(angio_fname);
angio = angio.I_seg_masked;
volshow(angio);

%% Smooth + Process segmentation prior to skeletonizing
%%% Processing Parameters
% Gaussian filter sigma
gsigma = 0.7;
% Minimum number of connected voxels
N = 50;

%%% Gaussian filter
angio_gaussian = imgaussfilt(angio, gsigma, 'FilterDomain', 'spatial');
% volshow(angio_gaussian);

%%% remove connected objects w/ fewer than N voxels
% Find connected components
CC = bwconncomp(angio_gaussian);
% Find number of voxels for each object
nvox = cellfun(@numel,CC.PixelIdxList);
% Find objects with fewer than N voxels
[~,idx] = find(nvox < N);
% Remove these objects
for ii=1:length(idx)
    angio_gaussian(CC.PixelIdxList{idx(ii)}) = 0;
end
% Show object after smoothing and removing
volshow(logical(angio_gaussian))
segmat2tif(uint8(angio_gaussian), fullfile(fout, 'angio_post.tif'))

%% Remove components that approximate spheres
% Use region properties to find objects with similar primary/secondary axes

%% Find region properties
% CC = bwconncomp(angio_gaussian);
% stats = regionprops3(angio_gaussian,'PrincipalAxisLength');


%% Skeletonize -- find branch points, end points, connected components
% TODO: find length of each segment within each component
% - subtract branch points from volume
% - perform bwconncomp
% - iterate through list of branch points (from same CC)
%   - add single branch point
%   - determine if new structure is connected

% Skeletonized unfiltered angio
angio_skel_pre = bwskel(logical(angio), "MinBranchLength",10);
volshow(angio_skel_pre);
segmat2tif(uint8(angio_skel_pre), fullfile(fout, 'angio_pre_skel.tif'))

% Convert 3D vol to skeleton
angio_skel_post = bwskel(logical(angio_gaussian), "MinBranchLength",10);
volshow(angio_skel_post);
segmat2tif(uint16(angio_skel_post), fullfile(fout, 'angio_post_skel.tif'))

%%% find branch points, end points, connected components
cc_skel = bwconncomp(angio_skel_post);
% retrieve length of each component
comlen = cc_skel.PixelIdxList;

% find branch points
skel_bp = bwmorph3(angio_skel_post, "branchpoints");

% find end points
skel_ep = bwmorph3(angio_skel_post, "endpoints");

% find connected components
label = bwlabeln(angio_skel_post);


%% Visualize loops for debugging
%{
volshow(angio_skel_post(200:300,50:125,1:25));
volshow(angio_gaussian(200:300,50:125,1:25));

volshow(angio_skel_pre(70:130,1:100,1:25));
volshow(angio(70:130,1:100,1:25));
%}
%% Find/remove loops from skeleton
% - iterate through list of connected components
%   - create volume for each separate connected component
%       - determine if contains end points
%       - if no end points, then it's a loop

% Iterate through list of connected components
for ii = 1:cc_skel.NumObjects   
    % 3D matrix of zeros, same size as original volume. This will be used
    % to visualize each connected component separately.
    vol_cc = zeros(cc_skel.ImageSize);
    
    % Convert only voxels in connected component to 1
    idx = cc_skel.PixelIdxList{ii};
    vol_cc(idx) = 1;
    
    %%% Determine if segment has end points
    seg_ep = bwmorph3(vol_cc, "endpoints");
    % If no endpoints, set component equal to zero because it's a loop
    if max(seg_ep(:))==0
        sprintf('Loop in connected component %i', ii)
        volshow(vol_cc)
        % show volume before removing
        % volshow(angio_skel)
        
        % Set components equal to zero 
        angio_skel_post(idx) = 0;
        
        % show volume after removing
        % volshow(angio_skel)  
    end
    
end

volshow(angio_skel_post)

%% Skeletonize and convert to Graph
% This section calls the function at the bottom of the script. This
% function then calls another function "seg_to_graph" which actually
% performs the skeletonization and converts the skeleton to the three
% dimensional graph (nodes + edges).
%{
segment_graph_data = seg_graph_init(segment, vox_dim, fullpath, filename);
%}
%% Calculate the vessel metrics from the graph data
% Metrics: length of each vessel, length density of volume , tortuosity
% NOTE: This section may require modification, depending on how you save
% your data in the previous section. 
%{
% Load graph and segmentation (angio)
graph = segment_graph_data.Graph; % might be segment_graph_data.Data.Graph
seg = segment_graph_data.angio;

% Load end nodes from graph
endnodes = graph.segInfo.segEndNodes;

% Load position of all nodes from graph
nodepos = graph.segInfo.segPos;

%%% Calculate total length (microns) from graph
len = graph.segInfo.segLen_um;
len_tot = sum(len(:));
met.total_length = len_tot;

%%% Calculate mean length (microns)
len_avg = mean(len);
met.avg_length = len_avg;

%%% Calculate length density
len_density = len_tot ./ vol;
met.length_density = len_density;

%%% Total number of vessels
nves = length(len);
met.total_vessels = nves;

%%% tortuosity arc-chord ratio (curve length / euclidean)
% Initalize matrix for storing tortuosity
tort = zeros(nves, 1);
for j=1:nves
    % convert segment end nodes to cartesian coordinate
    node1 = graph.nodes(endnodes(j,1), :);
    node2 = graph.nodes(endnodes(j,2), :);
    % Calcualte euclidean distance of segment
    euc = sqrt((node1(1) - node2(1)).^2 +...
                (node1(2) - node2(2)).^2 +...
                (node1(3) - node2(3)).^2);
    % Calculate tortuosity (single arc-chord ratio)
    tort(j) = len(j) ./ euc;
end
% Add to metrics structures
met.tortuosity = mean(tort);   

% Save the output
save(ad_cte_fout, 'met', '-v7.3');
%}
%% Function to skeletonize and convert to Graph
function [Data] = seg_graph_init(seg, vox_dim, fullpath, fname_seg)
% Initialize graph from segmentation
% INPUTS:
%   seg (mat): segmentation matrix
%   vox_dim (array): 3-element array of voxel dimensions (microns)
% OUTPUTS:
%   Data (struct): graph data for the segmentation volume

%%% Convert segmentation to graph (just the nodes and segments)
graph_nodes_segs = seg_to_graph(seg, vox_dim);

%%% Initialize graph metadata (Graph.Data)
[Data] = init_graph(graph_nodes_segs);

%%% Append "angio" data (segmentation matrix)
% Rearrange [x,y,z] --> [z,x,y]. This is the standard in
% the graph validation GUI.
angio = permute(seg, [3,1,2]);
Data.angio = angio;

% Create new filename for graph and add .MAT extension
fname_graph = strcat(fname_seg, '_graph_data.mat');
fout = fullfile(fullpath, fname_graph);
save(fout,'Data', '-v7.3');
end