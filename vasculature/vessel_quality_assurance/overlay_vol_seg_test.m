%% Main test script to overlay segmentation and original volume.
%{
The purpose of this script is to overlay the vessel segmentation and the
PS-OCT volume (each slice).
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

%% Initialize data path for linux or personal machine
% Check if running on local machine for debugging or on SCC for processing
if ispc
    % Top-level directory
    dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\';
    
    % Volume directory + volume filename (same for each subject)
    voldir = '\dist_corrected\volume\';
    volname = 'ref_4ds_norm_inv_crop2';

    % Directory to segmentation
    segdir = '\dist_corrected\volume\gsigma_1-3-5_gsize_5-13-21\';
    % Subject IDs
    subid = {'NC_6839'};
    
    % Filename to parse (this is test data)
    segname = 'ref_4ds_norm_inv_crop2_segment_pmin_0.23_mask40';

    % filename extension
    ext = '.tif';

elseif isunix
    error('This script was written to run just on windows.')
end

%% Overlay
%%% Load PSOCT volume (TIF) and convert to MAT
for ii = 1:length(subid)
    %%% Load segmentation
    % Define entire filepath
    fullpath = fullfile(dpath, subid{ii}, segdir);
    % Define filename of original ps-oct volume
    filename = strcat(fullpath, strcat(segname, ext));
    % Convert .tif to .MAT
    seg = TIFF2MAT(filename);
    
    %%% Load PS-OCT volume (inverted)
    fullpath = fullfile(dpath, subid{ii}, voldir);
    % Define filename of original ps-oct volume
    filename = strcat(fullpath, strcat(volname, ext));
    % Convert .tif to .MAT
    voli = TIFF2MAT(filename);
    
    %%% Make filepath
    % Save output volume
    fullpath = fullfile(dpath, subid{ii}, segdir);
    % Define filename for output of overlay
    filename = strcat(fullpath, strcat(segname, '_overlay', ext));
    
    %%% Call overlay function and save
    overlay_vol_seg(voli, seg, 'green', filename);
    
    %% XZ stack projection (for Stephan's paper)
%{
    %%% Range for cropping y-axis
    yrange = 200:300;

    %%% Minimum intensity projection of volume (XZ from non-inverted)
    % Convert from inverted to non-inverted
    vol = imcomplement(voli);
    % Crop Y to improve readability
    vol = vol(:,yrange, :);
    % Perform minimum intensity projection
    minp = min(vol,[],2);
    % Reshape into XZ
    [x,~,z] = size(minp);
    minp_re = reshape(minp, [x,z]);
    % Transpose matrix so that X coordinate is on x-axis of figure
    minp_re = minp_re';
    % Create figure for debugging
    figure; imagesc(minp_re); colormap(gray);
    set(gca,'YDir','normal');


    %%% Maximum intensity projection of segmentation (green)
    % Crop Y to improve readability
    seg = seg(:,yrange, :);
    % Perform minimum intensity projection
    maxp = max(seg,[],2);
    % Reshape into XZ
    [x,~,z] = size(maxp);
    maxp_re = reshape(maxp, [x,z]);
    % Transpose matrix so that X coordinate is on x-axis of figure
    maxp_re = maxp_re';
    % Create figure for debugging
    figure; imagesc(maxp_re); colormap(gray);
    set(gca,'YDir','normal');

    %%% Overlay projections
    ov = imoverlay(mat2gray(minp_re), maxp_re, 'green');
    figure; imshow(ov);
    % Update filename
    filename = strcat(filename(1:end-12),'_XZ_overlay.tif');
    imwrite(ov, filename);
%}
end
































