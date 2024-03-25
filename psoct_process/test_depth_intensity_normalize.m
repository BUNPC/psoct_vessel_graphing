%% Test script for the depth intensity normalization function
% TODO:
%{
- segment the WM in each slice
- create overlay of mask and tissue
- segment agarose and GM
- normalize across the z dimension for each layer of each physical slice.
- smooth the boundary between slices
%}

clear; clc; close all;


%% Add top-level directory of code repository to path
% This allows Matlab to find the functions in the project folders
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

%% Initialize data paths
% Path to top-level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
% Subfolder containing data
subdir = '/dist_corrected/volume/ref/';
% Filename of OCT volume (this will be the same for each subject)
volname = 'ref.tif';
% Filename of OCT intensity reference
irefname = 'ref_agarose_crop.tif';
% Subject ID list
subids = {'AD_10382', 'AD_20832', 'AD_20969',...
             'AD_21354', 'AD_21424',...
             'CTE_6489','CTE_6912',...
             'CTE_7019','CTE_8572','CTE_7126',...
             'NC_6047', 'NC_6839',...
             'NC_6974', 'NC_7597',...
             'NC_8095', 'NC_8653',...
             'NC_21499','NC_301181'};
subid = subids{3};

%% Import the OCT volume and intensity reference
%{
% Import OCT volume
subid = subids{3};
vol = fullfile(dpath, subid, subdir, volname);
vol = TIFF2MAT(vol);
% Import intensity reference volume
iref = fullfile(dpath, subid, subdir, irefname);
iref = TIFF2MAT(iref);
%}
%% Call normalization function
%{
% Normalize volume
volnorm = depth_intensity_normalize(vol, iref);
% Save output as a TIFF
fout = fullfile(dpath, subid, subdir, 'ref_agarose_norm.tif');
segmat2tif(volnorm,fout)
%}

%% Load the three normalized volumes
% Load the original non-normalized volume
ref = fullfile(dpath, subid, subdir, 'ref.tif');
ref = TIFF2MAT(ref);

% Load the volume normalized by white matter 
wm = fullfile(dpath, subid, subdir, 'ref_wm_norm.tif');
wm = TIFF2MAT(wm);

% Load the volume normalized by gray matter 
gm = fullfile(dpath, subid, subdir, 'ref_gm_norm.tif');
gm = TIFF2MAT(gm);

% Load the volume normalized by agarose 
ag = fullfile(dpath, subid, subdir, 'ref_agarose_norm.tif');
ag = TIFF2MAT(ag);

%% Measure variance across z stack

% Initialize struct for storing variance
v = struct();
% variable for number of slices above and below
n = 10;
% Store the wm, gm, and ag in struct
v(1).stack = wm;
v(2).stack = gm;
v(3).stack = ag;

for ii = 1:3
    % Create local variable of stack
    varmat = v(ii).stack;
    % Matrix for storing variance
    sub_var = zeros(size(varmat));    
    % Separate the z-stack into substacks of n frames
    for j = 1:size(wm,3)
        % Indices for slicing the variance matrix
        i_0 = max([1, j-n]);
        i_end = min([size(wm,3), j+n]);
        % Create subset of z-frames
        sub = varmat(:,:, i_0 : i_end);
        sub_var(:,:,j) = var(im2double(sub),0,3);
    end
    % Store variance in structure
    v(ii).var = sub_var;
end

%% Classify voxel by min variance
%%% Concatenate z-axis into two dimensional matrix
% Create local variables for white matter, gray matter, and agarose
wm_v = v(1).var;
gm_v = v(2).var;
ag_v = v(3).var;
% Reshape each into a 2d matrix
wm_r = reshape(wm_v, [size(wm_v,1), size(wm_v,2) .* size(wm_v,3)]);
gm_r = reshape(gm_v, [size(gm_v,1), size(gm_v,2) .* size(gm_v,3)]);
ag_r = reshape(ag_v, [size(ag_v,1), size(ag_v,2) .* size(ag_v,3)]);
% Stack all three
varstack = cat(3,wm_r, gm_r, ag_r);

%%% Compare variance across matrices (z-dimension)
% Find minimum across z-axis and assign to mask
[~,map] = min(varstack, [], 3);
% Reshape the mapping matrix back into the z-stack
map_r = uint8(reshape(map, [size(wm_v,1), size(wm_v,2), size(wm_v,3)]));
% Save the mapping
save(fullfile(dpath, subid, subdir, 'mask_map.mat'),'map_r');
%% Segment the background
% The edges are sometimes classified as white matter rather than agarose.
% They should be reclassified as agarose.
% TODO:
%   - serialize this ang run on each frame
%   - save mask output as TIF
%   - segment the agarose

%%% Load the mask mapping matrix
% load(fullfile(dpath, subid, subdir, 'mask_map.mat'));
% Take single slice for debugging 
bw = map_r == 1;
slice_n = 263;
bw = bw(:,:,slice_n);
figure; imagesc(bw); title('Before Processing')

%%% Detect large groups of pixels
% Find connected groups of pixels equal to 1
cc = bwconncomp(bw);
% Find k largest groups of pixels equal to 1
n_pix = cellfun(@numel, cc.PixelIdxList);
k = 4;
[~, idx] = maxk(n_pix,k);

%%% Find the group indices for pixels along the border.
% Subscripts of corners
c1 = [1, 1];
c2 = [size(bw,1), 1];
c3 = [1, size(bw,2)];
c4 = [size(bw,1), size(bw,2)];
% Convert corner subscripts to indices
cidx = ones(1,4);
cidx(1) = sub2ind(size(bw),c1(1), c1(2));
cidx(2) = sub2ind(size(bw),c2(1), c2(2));
cidx(3) = sub2ind(size(bw),c3(1), c3(2));
cidx(4) = sub2ind(size(bw),c4(1), c4(2));
% Array to track the pixel indices around the corners
border_idcs = [];
% Find groups of pixels along the border
for ii = 1:k
    % Extract pixel indices in group k
    p_idx = cc.PixelIdxList{idx(ii)};
    % Compare pixel indices to corner indices
    if any(ismember(cidx, p_idx))
        border_idcs = [border_idcs; p_idx];
    end        
end
% Convert boundary region to zeroes
bw(border_idcs) = 0;
figure; imagesc(bw); title('Border Removed');

%% Retain just the white matter
% % Dilate image
% se = strel('disk',1);
% bw = imdilate(bw,se);
% figure; imagesc(bw); title('Dilated');

% Find connected groups of pixels equal to 1
cc = bwconncomp(bw);
% Find largest connected group of pixels == 1
n_pix = cellfun(@numel, cc.PixelIdxList);
[~,idx] = max(n_pix);
largest_idcs = cc.PixelIdxList{idx};
% Discard all other pixelsw
del = 1:(size(bw,1) * size(bw,2));
del(largest_idcs) = [];
bw(del) = 0;
figure; imagesc(bw); title('Largest Group Retained')

% Fill image
bw = imfill(bw,'holes');
figure; imagesc(bw); title('Filled')

% Erode image to undo the dilation
se = strel('disk',1);
bw = imerode(bw,se);
figure; imagesc(bw); title('Eroded');

%%% Remove small groups of pixels
% cc = bwconncomp(bw);
% n_pix = cellfun(@numel, cc.PixelIdxList);
% p_idx = cc.PixelIdxList;
% del = p_idx(n_pix < 200);
% del = vertcat(del{:});
% bw(del) = 0;
% figure; imagesc(bw); title('small groups removed')
%% Segment the GM
% bwconncomp: find the largest WM region (1)

% Expand the WM region

% Imfill the WM region

% Remove the smaller sections of WM

%% Visual verification
% Plot the mapping matrix where: 1 = WM. 2 = GM. 3 = agarose
for ii = 1:20:264
    % Take slice from mapping matrix
    frame = map_r(:,:,ii);
    % Plot OCT frame and respective map
    figure;
    subplot(1,2,2); imshow(ref(:,:,ii)); colorbar
    subplot(1,2,1); imagesc(frame); colorbar
    % Title
    tstring = strcat('Frame ',num2str(ii));
    title(tstring);
end














