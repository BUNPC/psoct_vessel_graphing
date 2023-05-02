%% script for segmenting vessel .TIF
% Original Author: Jiarui Yang
% 10/21/20
% Modified by Mack Hyman (3/15/2023)
% Description:
%{
This is meant to be run on a local machine rather than the SCC. This is to
expedite debugging.

This takes the inverted image file as an input, applies a Gaussian filter,
then applies the Frangi filter. The user must select the parameters "thres"
and "sigma" in the section "Segment the inverted image."
%}

% To Do:
%{
- reapply mask to frangi filtered image
- turn this file into a function. then create a main script that calls this
        function.
%}

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

%% load inverted volume
fname = 'volume_ori_inv_cropped';

% Laptop directory structure
laptop_path = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_human_brain_resources\test_data\Hui_Frangi_dataset\200218depthnorm\';
filename = strcat(laptop_path, strcat(fname,'.tif'));
% Convert .tif to .MAT
ref = TIFF2MAT(filename);

%%% Quality assurance, view XZ profile
% cropped portion of volume
% view3D(ref)
% % Complete volume
% fname = 'ref_inv.tif';
% filename = strcat(laptop_path, fname);
% % Crop to be symmetric
% ref_complete = TIFF2MAT(filename);
% ref_complete = ref_complete(1:682,1:682,:);
% view3D(ref_complete)


%% Segment the inverted image
% I = 3D angiogram
% sigma = vector of standard deviation values of gaussian filter to
% thres = threshold value to create mask. Voxels are in grayscale [0, 1].
%           voxel < thres -->  0 (tissue)
%           voxel >= thres --> 1 (vessel)
%       "thres" is the limit to binarize (0=black, 1= white)

% Set threshold for cutoff
thres = 0.4;
% Set standard deviation for Gaussian filter kernel.
sigma = 2;

% Convert inverted image to double (required for segmentation).
I = double(ref);
[~, I_seg] = vesSegment(I, sigma, thres);
view3D(I_seg);

%% To Do: Apply mask to the segmented volume.
% The Gaussian filter may spread signal beyond the original masked area.
% This would result in an incorrectly labeled voxel in the binary matrix.
% I_seg_masked = 

%% Quality Assurance: view XY, XZ slices
%{
figure;
% X-Y slice original
subplot(2,2,1); imshow(ref(:,:,100)); title({'X-Y original','(slice 100)'})
% X-Y slice segment
subplot(2,2,3); imshow(I_seg(:,:,100)); title({'X-Y segment','(slice 100)'})
% X-Z slice original
subplot(2,2,2); imshow(squeeze(ref(:,100,:))); title({'X-Z original','(slice 100)'})
% X-Z slice segment
subplot(2,2,4); imshow(squeeze(I_seg(:,100,:))); title({'X-Z segment','(slice 100)'})
 text( 0.5, 1, fname, 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top' ) ;
%}

%% Save Vessel Segmentation
% Save vessel segment stack as .MAT for the next step (graph recon)
fout = strcat(laptop_path, fname, '_sigma', num2str(sigma));
save(strcat(fout,'.mat'), 'I_seg', '-v7.3');

% Save as .TIF for visual validation
segmat2tif(I_seg, strcat(fout,'.tif'));

