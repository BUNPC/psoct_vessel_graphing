%% Main script for segmenting 2PM images
% Author: Mack Hyman
% Email: mhyman@bu.edu
% Date Created: May 15, 2023
%
% Detailed Description
%{
This script performs the following:
- segment the original volume
- apply a mask to the segmentation
%}
clear; clc; close all;

%% Threshold the 2PM segmentations (debugging)
%{
minp = [0.2, 0.4, 0.6, 0.8, 0.9];
dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\AD_10382_2P\aip';
file = 'seg_sigma_1_3';
seg = load(fullfile(dpath, strcat(file,'.mat')));
seg = seg.seg;

for i=1:length(minp)
    % copy the original 2D frangi filter output
    tmp = seg;
    % Set anything below threshold equal to zero
    tmp( seg < minp(i)) = 0;
    % Save output as .TIF
    fname = strcat(file, 'minp_', num2str(minp(i)));
    fout = fullfile(dpath, strcat(fname, '.tif'));
    segmat2tif(tmp, fout);
end
%}
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

%% Import volume (.TIF or .BTF) & convert to MAT 

% Check if running on local machine for debugging or on SCC for processing
if ispc
    %%% Local machine
    dpath = ['C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\'...
        'test_data\BA4445_samples_16T'];
    % Subject IDs
    subid = {'\BA4445_I38_2PM'};
    subdir = '\aip\';
    % Number of slices to parse
    slices = [10];
    % Filename to parse (this is test data)
    fname = 'channel1';
    % filename extension
    ext = '.tif';
elseif isunix
    %%% Computing cluster (SCC)
    % Path to top-level directory
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_10T/';
    % Subject IDs
%     subid = {'AD_10382_2P', 'AD_20832_2P', 'AD_20969_2P', 'AD_21354_2P', 'AD_21424_2P',...
%              'CTE_6489_2P', 'CTE_6912_2P', 'CTE_7019_2P', 'CTE_8572_2P', 'CTE_7126_2P',...
%              'NC_21499_2P', 'NC_6047_2P', 'NC_6839_2P', 'NC_6974_2P', 'NC_7597_2P',...
%              'NC_8095_2P', 'NC_8653_2P'};
    subid = {'AD_10382_2P', 'AD_20832_2P'};
    subdir = '/aip/';
    % Filename to parse (this will be the same for each subject)
    fname = 'channel1';
    % filename extension
    ext = '.tif';
end

%% Frangi 2P Segmentation (Mathworks fibermetric)

% vessel thickness (pixels)
vth = 8;

% Minimum number of voxels to classify as segment
% A segment with < "min_conn" voxels will be removed
min_conn = 30;

% Length of line in the structuring element (strel) for expanding lines
line_len = 3;

% Iterate subjects
for ii = 1:length(subid)
    % Iterate slices
    for j = 1:length(slices)
        % Define filename for slice
        slice_fname = strcat(fname,'-',num2str(slices(j)),'_16b_cropped');
        % Define entire filepath 
        subpath = fullfile(dpath, subid{ii}, subdir);
        filename = strcat(subpath, strcat(slice_fname, ext));
        % Make sub-directory for slice
        slicepath = fullfile(subpath, strcat('slice_', num2str(slices(j))));
        if not(isfolder(slicepath))
            mkdir(slicepath)
        end
        % Convert .tif to .MAT
        vol = TIFF2MAT(filename);
        % Perform segmentation
        volmask = fibermetric(vol);
        % Dilate vessel mask
        [volmask_clean, fout] = clean_dilate(volmask, min_conn, line_len, subpath, slice_fname);
    end
end

%% Plot output of fibermetric
% Plot frangi output
figure; imshow(volmask)

% Binarize and overlay mask
bw = imbinarize(volmask);
figure; imshow(bw);
figure; imshow(labeloverlay(vol, bw));

%% Frangi 2P Segmentation (custom scripts)
%%% 2P microscopy voxel will always be [5, 5] micron
vox_dim = [2,2];

%%% Std. Dev. for gaussian filter (range)
% The value of sigma corresponds to the smallest resolvable radius of the
% vessels. If sigma equals 1, then the smallest resolvable vessel will be
% 1*voxel. For 2PM, the smallest resolvable vessel has a radius = 5um
sigma = [1, 3];

%%% Minimum number of voxels to classify as segment
% A segment with < "min_conn" voxels will be removed
min_conn = 30;

%%% Lenght of line in the structuring element (strel) for expanding lines
line_len = 3;

%%% Boolean for converting segment to graph (0 = don't convert, 1 = convert)
graph_boolean = 0;

%%% Segment each slice for each subject
for ii = 1:length(subid)
    for j = 1:length(slices)
        % Define filename for slice
        slice_fname = strcat(fname,'-',num2str(slices(j)),'_cropped_inv');
        % Define entire filepath 
        subpath = fullfile(dpath, subid{ii}, subdir);
        filename = strcat(subpath, strcat(slice_fname, ext));

        % Make sub-directory for slice
        slicepath = fullfile(subpath, strcat('slice_', num2str(slices(j))));
        if not(isfolder(slicepath))
            mkdir(slicepath)
        end
        % Convert .tif to .MAT
        vol = TIFF2MAT(filename);
        % Perform segmentation
        [seg, fname_seg] = segment_main(vol, sigma, slicepath, slice_fname);

        %%% Filter output of Frangi filter
        for k=1:length(line_len)
            % Remove small segments and apply mask
            [seg_masked, fout] = ...
                clean_dilate(seg, min_conn, line_len, slicepath, slice_fname);
            %%% Convert masked segmentation to graph
            if graph_boolean
                % Use masked segmentation to create graph
                Graph = seg_to_graph(seg_masked, vox_dim);
                
                % Initialize graph information (work in progress)
                % Graph = initialize_graph(Graph);
        
                % Create new filename for graph and add .MAT extension
                fname_graph = strcat(fname_seg,'_mask', num2str(radii(k)),'_graph.mat');
                fout = strcat(slicepath, fname_graph);
                save(fout,'Graph');
            end
        end
    end
end

%% Initialization of vesGraphValidate
function [graph_init] = initialize_graph(Graph)
%%% Perform the manual operations for initializing data in the GUI.
% Run "Verification > get segment info > Update"
% Run "Update branch info"
% Run "Regraph Nodes" to down sample
% Open GUI with both image and data (graph)
% Run prune_loops and prune_segment
% Run straighten
graph_init = Graph;
end

%% Apply Mask
function [seg_cleaned, fout] =...
    clean_dilate(seg, min_conn, line_len, fullpath, fname)
% Remove the edges labeled as vessels.
%   INPUTS:
%       seg (matrix) - output of segmentation function
%       min_conn (double) - minimum number of connections to classify as
%                           vessel
%       line_len (double) - length of line for strel
%       fullpath (string) - absolute directory for saving processed data
%       fname (string) - filename prior to applying mask
%   OUTPUTS:
%       seg_masked (matrix) - seg with boundaries eroded to remove
%           erroneously labeled vessels.

%%% Remove segments with fewer than "min_conn" connected voxels
figure; imshow(seg); title('Frangi - before processing')
[seg, ~, ~] = rm_small_seg(seg, min_conn);
figure; imshow(seg); title('Frangi - small segments removed')

%%% Dilate the lines
% Create perpendicular structure
se90 = strel('line',line_len,90);
se0 = strel('line',line_len,0);
% Dilate gradient mask
seg_dil = imdilate(seg, [se90 se0]);
figure; imshow(seg_dil); title('Frangi - Dilated Gradient Mask')

%%% Fill interior gaps
seg_fill = imfill(seg_dil,'holes');
figure; imshow(seg_fill); title('Frangi - Filled Holes')

%%% Smooth segmentation
seD = strel('diamond',1);
seg_final = imerode(seg_fill,seD);
seg_final = imerode(seg_final,seD);
figure; imshow(seg_final); title('Frangi - final');

%%% Save segmented/masked volume as .MAT and .TIF
% Convert processed matrix to tif
tmp_fname = strcat(fname,'_processed');
fout = strcat(fullpath, tmp_fname, '.tif');
segmat2tif(seg_cleaned, fout);

% Save vessel segment stack as .MAT for the next step (graphing)
fout = strcat(fullpath, tmp_fname, '.mat');
save(fout, 'seg_cleaned', '-v7.3');

end

%% Segment volume
function [seg, fname] =...
    segment_main(vol, sigma, fullpath, fname)
% Multiscale vessel segmentation
%   INPUTS:
%       vol (matrix) - the original volume prior to segmentation
%       FrangiScaleRange (array) - range of std. dev. values of gaussian
%            filter to calcualte hessian matrix at each voxel
%       min_prob - threshold to determine which voxel belongs to a vessel.
%           Applied to probability matrix from frangi filter output
%       min_conn - vesSegment uses the function bwconncomp to determine the
%               number of connected voxels for each segment. If the number
%               of voxels is less than this threshold, then the segment
%               will be removed.
%       fullpath (string) - absolute directory for saving processed data
%       fname (string) - filename of segmentation
%       ext (string) - filename extension (.TIF or .MAT)
%   OUTPUTS:
%       seg - segmentation of vessels (Frangi filter of original volume)
%       fname (string) - upadted filename prior to applying mask

%% convert volume to double matrix
I = double(vol);

%%% Range of sigma values
options.FrangiScaleRange = sigma;

%%% Segment the original volume
[seg, ~, ~] = FrangiFilter2D(I, options);

%%% Save segmentation
fname = strcat(fname,'_segment','_sigma_',...
    num2str(sigma(1)), '_', num2str(sigma(2)));
fout = fullfile(fullpath, strcat(fname, '.tif'));
segmat2tif(seg, fout);

end