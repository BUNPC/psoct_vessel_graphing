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

clear; clc;

%% Imaging parameters
% voxel dimensions
vox_dim = [10, 10, 10]; % what are the units here? Leftover from Jairui.

%% Segment volume 
% Load masked volume
fname = 'volume_nor_inverted_masked';
% Laptop directory structure
laptop_path = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_human_brain_resources\test_data\Hui_Frangi_dataset\200726PSOCT\';
filename = strcat(laptop_path, strcat(fname,'.tif'));
% Convert .tif to .MAT
vol = TIFF2MAT(filename);

% Call segmentation
Ithresh = 0.2;

%% Convert segmentation to graph

%% Initialization of vesGraphValidate

%% Perform manual pruning
% The user must use the matlab GUI to manually remove these segments

%% Calculate Diameter
% Load Graph struct
fname = 'volume_nor_inverted_masked_sigma1_frangi_seg_regraphed';
laptop_path =...
    'C:\Users\mack\Documents\BU\Boas_Lab\psoct_human_brain_resources\test_data\Hui_Frangi_dataset\200726PSOCT\';
fpath = strcat(laptop_path, strcat(fname,'.mat'));
Data = load(fpath,'Graph');

% Call function to calculate diameter at each node
Diam = GetDiam_graph(...
    vol,...
    Data.Graph.nodes,...
    Data.Graph.edges,...
    Ithresh,...
    vox_dim);

%% Calculate Tortuosity
tortuosity = vessel_tortuosity_index(Data.Graph, Ithresh);

%% Histograms for geometries
% Histo for diameter
figure; histogram(Diam);
title('Vessel Diameter')
xlabel('Diameter (microns)')
ylabel('Count')
set(gca, 'FontSize', 20)

% Histo for diameter
figure; histogram(tortuosity);
title('Vessel Tortuosity')
xlabel('Tortuosity (unitless)')
ylabel('Count'); ylim([0,200])
set(gca, 'FontSize', 20)


%% (OLD CODE) Median diameter
%{
dia_vessel=zeros(1,length(Data.Graph.segInfo.segLen));
for i=1:length(dia_vessel)
    dia_vessel(i)=median(Diam(find(Data.Graph.segInfo.nodeSegN(:)==i)));
end
figure; histogram(dia_vessel,'BinWidth',10);
%}

%% (OLD CODE) length of each vessel
%{
% remove fake vessels
length_vessel( length_vessel(:)==0 ) = [];
figure;
histogram(length_vessel.*10,0:25:1000);
%}


%% (OLD CODE) extract vessels and clean up boundaries by creating a mask
%{
V_seg=TIFF2MAT('I_seg_Ann.tif');
mask=TIFF2MAT('I_mask.tif');
for i=1:size(img,3)
    mask_tmp=squeeze(mask(:,:,i));
    mask_tmp(mask_tmp(:)~=0)=1;
    mask_tmp = bwareaopen(logical(mask_tmp), 500);
    V_seg(:,:,i)=uint16(mask_tmp).*squeeze(V_seg(:,:,i));
    % figure;imshow(mask_tmp);
end
MAT2TIFF(V_seg,'I_seg_masked2.tif');
%}

%% (OLD CODE) calculate topology
%{
% skeletonization
I_seg=TIFF2MAT('I_seg_nor_mask.tif');
I_seg(I_seg(:)~=0)=1;
I_skel=bwskel(logical(I_seg));
MAT2TIFF(I_skel,'I_seg_skel.tif');
%}