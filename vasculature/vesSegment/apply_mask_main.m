%% Apply the mask.tif to the volume.tif
% In some cases, the original volume still contains the background signal.
% In the case that there exists a mask.tif file for this volume, this
% script will call a function to apply the mask to the volume.

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

%% Lookup table of pia pixel intensity reference values
subid = {'AD_10382', 'AD_20832', 'AD_20969',...
             'AD_21354', 'AD_21424',...
             'CTE_6489','CTE_6912',...
             'CTE_7019','CTE_8572','CTE_7126',...
             'NC_6047', 'NC_6839',...
             'NC_6974', 'NC_7597',...
             'NC_8095', 'NC_8653',...
             'NC_21499','NC_301181'};

subid =     {'NC_6047', 'NC_6839',...
             'NC_6974', 'NC_7597',...
             'NC_8095', 'NC_8653',...
             'NC_21499','NC_301181'};

subid = {'AD_20969','AD_21424'};

piaref = struct;
piaref.AD_10382 = 6e4;
piaref.AD_20832 = 6e4;
piaref.AD_20969 = 6e4;
piaref.AD_21354 = 6e4;
piaref.AD_21424 = 6e4;
piaref.CTE_6489 = 6e4;
piaref.CTE_6912 = 6e4;
piaref.CTE_7019 = 6e4;
piaref.CTE_8572 = 6e4;
piaref.CTE_7126 = 6e4;
piaref.NC_6047 = 6e4;
piaref.NC_6839 = 6e4;
piaref.NC_6974 = 6e4;
piaref.NC_7597 = 6e4;
piaref.NC_8095 = 6e4;
piaref.NC_8653 = 6e4;
piaref.NC_21499 = 6e4;
piaref.NC_301181 = 6e4;

%% Import files
if ispc
    % Laptop directory structure
    dpath = ['C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\' ...
            'test_data\Ann_Mckee_samples_10T\'];
    % Subfolder containing data
    subdir1 = '\dist_corrected\volume\';
elseif isunix
    % Upper level directory
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T';    
    % Subfolder containing data
    subdir1 = '/dist_corrected/volume/';
end

% Filename with PSOCT scattering tissue volume
vol_name = 'ref_4ds_norm_inv'; 
% Combined segmentations subfolder
subdir2 = 'combined_segs';
% Combined segmentation subfolder
segdir = 'gsigma_3-5-7_5-7-9_7-9-11';
% Filename with combined segmentations
seg_name = 'seg';


%% Iterate through subjects. Create and apply mask.
for ii = 1:length(subid)
    % Pixel intensity of pia
    ref = piaref.(subid{ii});

    % Import volume
    fpath = fullfile(dpath, subid{ii}, subdir1, strcat(vol_name,'.tif'));
    vol = TIFF2MAT(fpath);
    
    % Import combined segmentation file
    fpath = fullfile(dpath, subid{ii}, subdir1, subdir2, segdir, strcat(seg_name,'.mat'));
    seg = load(fpath,'seg');
    seg = seg.seg;
    
    % Call function to mask background and pia
    [volm, segm] = create_mask(vol, seg, ref);
       
    % Export masked volume and segmentation
    volm_out = fullfile(dpath, subid{ii}, subdir1, strcat(vol_name,'_masked.tif'));
    save(volm_out, "volm", '-v7.3');
    
    % Export masked segmentation
    segm_out = fullfile(dpath, subid{ii}, subdir1, subdir2, segdir, strcat(seg_name,'_masked.mat'));
    save(segm_out, "segm", '-v7.3');

    % Overlay the segmentation and masked volume (save overlay)
    fout = fullfile(dpath, subid{ii}, subdir1, subdir2, segdir, strcat(seg_name, '_masked_overlay.tif'));
    overlay_vol_seg(volm, segm, 'green', fout);
    
end

%% Debugging: View single slice for reference
%{
slice = vol(:,:,100);
figure; imshow(slice)
% Create boundary mask
mask = boundarymask(slice);

% Create structuring element
[m,n] = size(slice);
se = strel('diamond',20);
% Erode the mask
maskeroded = imerode(mask,se);

% Plot the mask and eroded mask
figure;
subplot(1,2,1); imshow(mask); title('boundary mask')
subplot(1,2,2); imshow(maskeroded); title('Eroded')

% Overlay mask and original
% slice_me = slice .* uint16(maskeroded);
% figure; imshow(slice_me); title('Mask eroded')

%}
%% Debugging: Find gray difference weight for removing pia
% The grraydiffweight can find pixels with a similar pixel intensity as the
% reference. By setting the reference equal to a pixel in the pia, this
% function will find all the pia along the border. However, it will also
% falsely label some vessel lumens as pia.
%{
% Find difference of gray values
ref = 6e4;
w = graydiffweight(slice, ref, 'GrayDifferenceCutoff', 10000);

% Threshold the difference map
wth = log(w);
wth(wth<4) = 0;
figure; imshow(wth);

% Dilate the segmentation
se = strel('disk',5);
wth = imdilate(wth, se);
figure; imshow(wth);

% Keep components with number of connections within range
wth = logical(wth);
range = [1e4, 1e8];
bw = bwareafilt(wth, range);
% Dilate again
bw = imdilate(bw, se);
figure; imshow(bw);

% Invert the mask
piamask = ~bw;
figure; imshow(piamask); title('Piamask before connectivity')
% Keep components with number of connections within range
piamask = bwareafilt(piamask, range);
figure; imshow(piamask); title('Piamask after connectivity')
%}
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
fout = strcat(dpath, strcat(vol_name,'_masked_eroded_island_rm.tif'));
segmat2tif(vol_masked, fout);
%}


