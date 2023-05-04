%% Main file for calling segmentation functions
% Author: Mack Hyman
% Date Created: March 16, 2023
%
% Detailed Description
%{
This script performs the following:
- segmentation ()
- convert segmentation to graph ()
- prune graph (remove loops and unterminated segments)
    - remove loops ()
    - remove segments ()
- overlay graph and image
- diameter
- tortuosity (vessel_tortuosity_index.m)
- length
- topology (why does this use the mask?)
%}
clear; clc; close all;

%% Add top-level directory of code repository to path
% This allows Matlab to find the functions in the project folders
% Print current working directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Remove the two sub folders to reach top-level directory
% (psoct_vessel_graphing)
topdir = mydir(1:idcs(end));
addpath(genpath(topdir));


%% Import volume (.TIF or .BTF) & convert to MAT 

%%% Local paths (windows PC)
%{
dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_human_brain_resources\test_data\Hui_Frangi_dataset\200218depthnorm\';
fname = 'volume_ori_inv_cropped';
% filename extension
ext = '.tif';
filename = strcat(dpath, strcat(fname,ext));
% Convert .tif to .MAT
vol = TIFF2MAT(filename);
%}

%%% SCC paths (windows PC)

dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_10T/AD_10382/dist_corrected/volume/';
fname = 'ref_4ds_norm';
% filename extension
ext = '.btf';
filename = strcat(dpath, strcat(fname,ext));
% Convert .tif to .MAT
vol = TIFF2MAT(filename);
%}

%% Set Voxel Size
%%% Assign PS-OCT voxel dimension [x, y, z] according to downsample factor
% Downasample factor = 4 --> Voxel = [12, 12, 15] micron
% Downasample factor = 10 --> Voxel = [30, 30, 35] micron
% TODO: Create function to extract downsample and voxel size from filename.
%       The filename contains either "_10ds" or "_4ds".

% Temporary hard-coded voxel dimensions
% PS-OCT voxel dimensions (microns)
vox_dim = [30, 30, 35];

%%% 2P microscopy voxel will always be 5um x 5um
vox2p = [5, 5];

%% Segment the volume
% Std. Dev. for gaussian filter
sigma = 1;
% threshold to determine whether voxel belongs to a vessel. This is applied
% to the probability matrix from the output of the frangi filter.
thresh = 0.2;
% Segments with few than "min_conn" voxels will be removed
min_conn = 30;

% Run segmentation function
[I_seg] = ...
    segment_main(dpath, fname, ext, sigma, thresh, min_conn);

% Save segmentation before applying mask
fname = strcat(fname,'_segment','_sigma', num2str(sigma));
fout = strcat(dpath, fname, '.tif');
segmat2tif(I_seg, fout);

% Call masking function

%% Apply mask to segmentation volume -- remove erroneous vessels
% TODO: find optimal range for remove_mask_islands
% TODO: create function "clean_mask" and perform:
%       - convert original volume to binary (backgroun = 0)
%       - imerode - shrink boundaries of mask
%       - remove_mask_islands (remove islands of pixels from mask)
%       - apply_mask (multiply mask and segmented volume)
%       - save output

%%% Create 3D mask
mask = logical(vol);

radii = 10:2:22;

for ii = 1:10
    %%% Erode mask to remove small pixels on border that are not part of volume
    se = strel('disk', radii(ii));
    mask = imerode(mask, se);
    
    %%% Remove islands of pixels from mask
    % Range of object size to keep
    range = [1e4, 1e8];
    mask = remove_mask_islands(mask, range);
    
    %%% Apply mask to segmentation volume
    % Convert from logical back to uint16 for matrix multiplication
    mask = uint16(mask);
    % Element-wise multiply mask and volume
    I_seg_masked = apply_mask(I_seg, mask);
    
    %%% Save segmented/masked volume as .MAT and .TIF
    % Convert masked image back to tif
    fname = strcat(fname,'_masked_radius_', num2str(radii(ii)));
    fout = strcat(dpath, fname, '.tif');
    segmat2tif(I_seg_masked, fout);
    % Save vessel segment stack as .MAT for the next step (graph recon)
    fout = strcat(dpath, fname, '.mat');
    save(fout, 'I_seg_masked', '-v7.3');
end

%% Convert segmentation to graph

% Use masked segmentation to create graph
Graph = seg_to_graph(fout, vox_dim);

% Create new filename for graph and add .MAT extension
fname = strcat(fname, '_graph.mat');
fout = strcat(dpath, fname);
save(fout,'Graph');

%% Initialization of vesGraphValidate
% Run "Verification > get segment info > Update"
% Run "Update branch info"
% Run "Regraph Nodes" to down sample
% Open GUI with both image and data (graph)
% Run prune_loops and prune_segment
% Run straighten

%% Apply Mask
function [I_seg_masked] = mask_segments(I, I_seg, epsilon)
% Remove the edges labeled as vessels.
%   INPUTS:
%       I (matrix) - original volume
%       I_seg (matrix) - output of segmentation function
%       epsilon (double) - scalar for determining masking
%   OUTPUTS:
%       I_seg_masked (matrix) - I_seg with boundaries eroded to remove
%           erroneously labeled vessels.

%%% Create mask from normalized volume
% TODO: convert volume -> binary 3D mask
mask = logical(I);

%%% Erode mask to remove small pixels on border that are not part of volume
se = strel('disk',10);
mask = imerode(mask, se);

%%% Remove islands of pixels from mask
% Range of object size to keep
range = [1e4, 1e8];
mask = remove_mask_islands(mask, range);

%%% Apply mask to segmentation volume
% Convert from logical back to uint16 for matrix multiplication
mask = uint16(mask);
% Element-wise multiply mask and volume
vol_masked = apply_mask(vol, mask);
% Convert masked image back to tif
fout = strcat(laptop_path, strcat(vol_name,'_masked.tif'));
segmat2tif(vol_masked, fout);
%}

end

%% Segment volume
function [I_seg] =...
    segment_main(dpath, fname, ext, sigma, thresh, min_conn)
% Multiscale vessel segmentation
%   INPUTS:
%   sigma - vector of standard deviation values of gaussian filter to
%           calcualte hessian matrix at each voxel
%   thres - threshold to determine which voxel belongs to a vessel. This is
%           applied to the probability matrix from the output of the frangi
%           filter.
%   min_conn - vesSegment uses the function bwconncomp to determine the
%               number of connected voxels for each segment. If the number
%               of voxels is less than this threshold, then the segment
%               will be removed.
%   OUTPUTS:
%       I_seg - segmentation of vessels (Frangi filter of original volume)
%       I_seg_masked - masked segmentation (removed borders)

%%% Load .TIF and convert to .MAT
filename = strcat(dpath, strcat(fname,ext));
vol = TIFF2MAT(filename);
% convert volume to double matrix
I = double(vol);

%%% Segment the original volume
[~, I_seg] = vesSegment(I, sigma, thresh, min_conn);

%%% Apply a mask to segmentation
% Scalar for determining how much to erode mask
% epsilon = 1;
% I_seg_masked = mask_segments(I, I_seg, epsilon);

end

%% Convert segmentation to graph
function graph_main()

end