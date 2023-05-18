%% Main file for calling segmentation functions
% Author: Mack Hyman
% Date Created: March 16, 2023
%
% Detailed Description
%{
This script performs the following:
- segment the original volume
- apply a mask to the segmentation
- convert segmentation to graph
To Do:
- find optimal range for remove_mask_islands
- prune graph (remove loops and unterminated segments)
    - remove loops ()
    - remove segments ()
%}
clear; clc; close all;
tic

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
    dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\';
    % Subject IDs
    subid = {'CTE_7019'};
    subdir = '\dist_corrected\volume\';
    % Filename to parse (this is test data)
    fname = 'ref_4ds_norm_inv_cropped';
    % filename extension
    ext = '.tif';
elseif isunix
    %%% Computing cluster (SCC)
    % Path to top-level directory
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_10T/';
    % Complete subject ID list for Ann_Mckee_samples_10T
%     subid = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424',...
%              'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_8572', 'CTE_7126',...
%              'NC_21499', 'NC_6047', 'NC_6839', 'NC_6974', 'NC_7597',...
%              'NC_8095', 'NC_8653'};
    % Partial subject ID list for testing script on SCC
    subid = {'AD_20969', 'AD_21354'};
    subdir = '/dist_corrected/volume/';
    % Filename to parse (this will be the same for each subject)
    fname = 'ref_4ds_norm_inv';
    % filename extension
    ext = '.btf';    
end

%% Initialization parameters

%%% Assign PS-OCT voxel dimension [x, y, z] according to downsample factor
% Downasample factor = 4 --> Voxel = [12, 12, 15] micron
% Downasample factor = 10 --> Voxel = [30, 30, 35] micron
% 2P microscopy voxel will always be 5um x 5um

% Set voxel dimensions from filename
if regexp(fname, '4ds')
    vox_dim = [12, 12, 15];
elseif regexp(fname, '10ds')
    vox_dim = [30, 30, 35];
else
    vox_dim = [30, 30, 35];
end

% Std. Dev. for gaussian filter (one value or array)
% The value of sigma corresponds to the smallest resolvable radius of the
% vessels. If sigma==1, then the smallest resolvable vessel will be
% 1*voxel. In our case, the smallest resolvable vessel has a radius = 12um
sigma = 1;

% Minimum fringi filter probability to classify voxel as vessel
min_prob = 0.20:0.02:0.26;
% A segment with < "min_conn" voxels will be removed
min_conn = 30;

% Array (or single value) of radii for eroding the mask
radii = 40;

% Boolean for converting segment to graph (0 = don't convert, 1 = convert)
graph_boolean = 0;

for ii = 1:length(subid)
    %% Segment the volume
    % Define entire filepath 
    fullpath = fullfile(dpath, subid{ii}, subdir);
    filename = strcat(fullpath, strcat(fname, ext));
    % Convert .tif to .MAT
    vol = TIFF2MAT(filename);
    
    for j = 1:length(min_prob)
        [I_seg, fname_seg] = ...
            segment_main(vol, sigma, min_prob(j), min_conn, fullpath, fname);
        
        %% Mask segmented volume (remove erroneous vessels) & Convert to Graph
        % The function for creating the mask requires a radius. This for-loop will
        % iterate over an array of radii. For each radius, it will create a mask,
        % apply the mask to the segmentation volume, and save the output.
        % If the graph_boolean is true (1), then the masked segmentation will be
        % converted to a graph.
    
        % Create 3D mask from original volume
        mask = logical(vol);
        for k = 1:length(radii)
            %%% Apply mask and save .MAT and .TIF
            [I_seg_masked] = mask_segments(I_seg, mask, radii(k), fullpath, fname_seg);
            
            %%% Convert masked segmentation to graph
            if graph_boolean
                % Use masked segmentation to create graph
                Graph = seg_to_graph(I_seg_masked, vox_dim);
                
                % Initialize graph information (work in progress)
                % Graph = initialize_graph(Graph);
        
                % Create new filename for graph and add .MAT extension
                fname_graph = strcat(fname_seg,'_mask', num2str(radii(k)),'_graph.mat');
                fout = strcat(fullpath, fname_graph);
                save(fout,'Graph');
            end
        end
    end
end
toc

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
function [I_seg_masked] = mask_segments(I_seg, mask, radius, fullpath, fname)
% Remove the edges labeled as vessels.
%   INPUTS:
%       I_seg (matrix) - output of segmentation function
%       mask (matrix) - unsegmented volume converted to logicals
%       radius (double array) - radius of disk for eroding the mask
%       fullpath (string) - absolute directory for saving processed data
%       fname (string) - filename prior to applying mask
%   OUTPUTS:
%       I_seg_masked (matrix) - I_seg with boundaries eroded to remove
%           erroneously labeled vessels.

%%% Erode mask to remove small pixels on border that are not part of volume
se = strel('disk', radius);
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
tmp_fname = strcat(fname,'_mask', num2str(radius));
fout = strcat(fullpath, tmp_fname, '.tif');
segmat2tif(I_seg_masked, fout);
% Save vessel segment stack as .MAT for the next step (graph recon)
fout = strcat(fullpath, tmp_fname, '.mat');
save(fout, 'I_seg_masked', '-v7.3');

end

%% Segment volume
function [I_seg, fname] =...
    segment_main(vol, sigma, min_prob, min_conn, fullpath, fname)
% Multiscale vessel segmentation
%   INPUTS:
%       vol (matrix) - the original volume prior to segmentation
%       sigma (array) - vector of std. dev. values of gaussian filter to
%           calcualte hessian matrix at each voxel
%       thres - threshold to determine which voxel belongs to a vessel.
%           Applied to probability matrix from frangi filter output
%       min_conn - vesSegment uses the function bwconncomp to determine the
%               number of connected voxels for each segment. If the number
%               of voxels is less than this threshold, then the segment
%               will be removed.
%       fullpath (string) - absolute directory for saving processed data
%       fname (string) - filename of segmentation
%       ext (string) - filename extension (.TIF or .MAT)
%   OUTPUTS:
%       I_seg - segmentation of vessels (Frangi filter of original volume)
%       fname (string) - upadted filename prior to applying mask

%% convert volume to double matrix
I = double(vol);

%%% Segment the original volume
[~, I_seg] = vesSegment(I, sigma, min_prob, min_conn);

%%% Save segmentation
fname = strcat(fname,'_segment','_sigma', num2str(sigma), '_thresh', num2str(min_prob));
fout = strcat(fullpath, fname, '.tif');
segmat2tif(I_seg, fout);

end