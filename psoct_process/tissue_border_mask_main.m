%% Create mask for tissue and pia boundary
% This script will find the edge of the tissue for each slice, create a
% continuous polygon to outline the tissue, and then erode the boundary of
% the tissue. This will be repeated for each slice. This mask will then be
% applied to the segmentation to remove false positives along the border of
% the tissue.

% TODO:
%{
- Add code to remove islands from mask
- make another create_mask function "create_mask_v0"
    - include threshold, active contour, island removal
- make this the main script, iterate over subjects, call:
    - "create_mask" for all subjects
    - "apply_mask" for all subjects
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

%%% Local machine
if ispc
    dpath = ['C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures'...
            '\test_data\Ann_Mckee_samples_10T\'];
%%% Computing cluster (SCC)
elseif isunix
    % Path to top-level directory
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
end

% Subfolder containing data
subdir = '/dist_corrected/volume/';
% Filename to parse (this will be the same for each subject)
fname = 'ref.tif';
%%% Complete subject ID list for Ann_Mckee_samples_10T
subid = {'AD_10382', 'AD_20832', 'AD_20969',...
         'AD_21354', 'AD_21424',...
         'CTE_6489','CTE_6912',...
         'CTE_7019','CTE_8572','CTE_7126',...
         'NC_6047', 'NC_6839',...
         'NC_6974', 'NC_7597',...
         'NC_8095', 'NC_8653',...
         'NC_21499','NC_301181'};

%% Load raw volume (TIF) and convert to MAT
% Define entire filepath 
fullpath = fullfile(dpath, subid{1}, subdir);
% Define file path
fpath = fullfile(fullpath, fname);
% Import the TIF stack
d = TIFF2MAT(fpath);
% Take slice
s = d(:,:,1); figure; imshow(s); title('Original');

%% Multithreshold
% Create multithreshold
lvl = multithresh(s,1);
maskth = imquantize(s, lvl);
figure; imshow(maskth, []); title('Multithresholded')
% Fill the image
mask_fill = imfill(maskth);
figure; imshow(maskth, []); title('Multithresholded & Filled')
%% Threshold
% Find quartile of median
imin = 0.25*median(s(:));
% Remove voxels below cutoff
sth = s;
sth(s<imin) = 0; figure; imshow(sth); title('Thresholded');
% Find edges
bw = edge(sth, "canny"); figure; imshow(bw); title('Thresholded & Edges');

%% Active contour
% Create starting mask
mask = zeros(size(s));
mask(50:1642,10:2600) = 1;

%%% Apply active contour
niter = 1e4;
bw = activecontour(sth, mask, niter);
figure; imshow(bw); title('Active Contour, 1e4 Iterations');

%%% Erode the border
se = strel('disk',10);
bwfill = imerode(bw, se);
figure; imshow(bwfill); title('Active Contour, 1e4 Iterations, Eroded 5');

%%% Remove islands of pixels


%%% Overlay tissue with active contour mask
% Binarize the mask ([0, 255] to [0, 1])
mask_bin = uint16(bwfill);
% Element-wise multiplication of mask and volume
masked = uint16(s) .* mask_bin;
% Convert masked from logical back to grayscale [0, 255]
masked(masked==1) = 255;
% Display masked tissue
figure; imshow(masked); title('Masked Tissue Slice')




