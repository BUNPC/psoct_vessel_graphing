%% Apply the mask.tif to the volume.tif
% In some cases, the original volume still contains the background signal.
% In the case that there exists a mask.tif file for this volume, this
% script will call a function to apply the mask to the volume.

%% TODO:
% - run with Hui_Frangi_dataset\200726PSOCT
% - convert TIF files to .MAT

clear; clc; close all;

%% Add top-level directory of code repository to path
% Print current working directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Remove the two sub folders to reach parent
% (psoct_human_brain\vasculature\vesSegment)
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));

%% Import files
if ispc
    % Laptop directory structure
    laptop_path = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_human_brain_resources\test_data\Hui_Frangi_dataset\200726PSOCT\';
    
    % Import volume
    vol_name = 'volume_nor_inverted';
    filename = strcat(laptop_path, strcat(vol_name,'.tif'));
    vol = TIFF2MAT(filename);
    
    % Import mask
    mask_name = 'volume_mask';
    filename = strcat(laptop_path, strcat(mask_name,'.tif'));
    mask = TIFF2MAT(filename);
    
    % Verify matrices of same size
    if size(vol) ~= size(mask)
        error('Mask and volume matrices are different sizes.')
    end

elseif isunix
    dpath = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/CTE_7126/' ...
        'dist_corrected/volume/gsigma_3-5-7_gsize_13-21-29/'];
    vol_name = 'ref_4ds_norm_inv_segment_pmin_0.22';
    filename = fullfile(dpath, strcat(vol_name,'.tif'));
    vol = TIFF2MAT(filename);
end

%% Segment


%% Perform operations on mask and image
%{
% TODO: find optimal range for remove_mask_islands
% TODO: create function "clean_mask" and perform both:
%       - imerode - remove boundaries
%       - remove_mask_islands - remove islands of pixels

%%% Erode mask to remove small pixels on border that are not part of volume
se = strel('disk',10);
mask = imerode(mask, se);

%%% Remove islands of pixels from mask
% Range of object size to keep
range = [1e4, 1e8];
mask_isl = remove_mask_islands(mask, range);

%%% Apply mask to volume
% Convert from logical back to uint16 for matrix multiplication
mask_isl = uint16(mask_isl);
% Element-wise multiply mask and volume
vol_masked = apply_mask(vol, mask_isl);
% Convert masked image back to tif
fout = strcat(laptop_path, strcat(vol_name,'_masked_eroded_island_rm.tif'));
segmat2tif(vol_masked, fout);
%}


