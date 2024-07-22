%% Main file for calling segmentation functions
% Author: Mack Hyman
% Date Created: June 18, 2024
%
% Detailed Description
%{
1) Test a range of bwskel minimum branch lengths
2) Create 3D overlays to visualize the results
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
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));

%%% Set maximum number of cores
if isunix
    % Retrieve the number of available cores
    n_cores = str2num(getenv('NSLOTS'));
    % Set the maximum number of threads equal to the number of cores
    maxNumCompThreads(n_cores);
end

%% Initialize data path for linux or personal machine
if isunix
    % Path to top-level directory
    dpath = '/projectnb/npbssmic/ns/CAA/';
    % Subfolder containing data
    subdir = 'segmentations/';
    % Filename to parse (this will be the same for each subject)
    seg_names = {'caa6-frontal_vessels-masked',...
              'caa6-occipital_vessels-masked',...
              'caa17_occipital_THRESH-0.5_masked',...
              'caa22-frontal_vessels-masked',...
              'caa22-occipital_vessels-masked',...
              'caa25-frontal_vessels-masked',...
              'caa25-occipital_vessels-masked',...
              'caa26-frontal_vessels-masked',...
              'caa26-occipital_vessels-masked'};
    % filename extension
    ext = '.mat';
    
    %%% Complete subject ID list for CAA samples
    subids = {'caa6/frontal/', 'caa6/occipital/',...
             'caa17/occipital/',...
             'caa22/frontal/','caa22/occipital/',...
             'caa25/frontal/', 'caa25/occipital/',...
             'caa26/frontal/', 'caa26/occipital/'};
elseif ispc
    % Path to top-level directory
    dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\Martinos_Datasets\';
    % Subfolder containing data
    subdir = 'segmentations\';
    % Complete subject ID list for CAA samples
    subids = {'caa_6/frontal/', 'caa_6/occipital/',...
             'caa_17/occipital/',...
             'caa_22/frontal/','caa_22/occipital/',...
             'caa_25/frontal/', 'caa_25/occipital/',...
             'caa_26/frontal/', 'caa_26/occipital/'};
    % Filename of vessel segmentation
    seg_names = {'caa6-frontal_vessels-masked.mat',...
              'caa6-occipital_vessels-masked.mat',...
              'caa17_occipital_THRESH-0.5_masked.mat',...
              'caa22-frontal_vessels-masked.mat',...
              'caa22-occipital_vessels-masked.mat',...
              'caa25-frontal_vessels-masked.mat',...
              'caa25-occipital_vessels-masked.mat',...
              'caa26-frontal_vessels-masked.mat',...
              'caa26-occipital_vessels-masked.mat'};
    % Name of the OCT volume
    mus_name = 'mus.nii';
    % Filename of EPVS segmentation
    epvs_dir = {'caa_17\occipital\segmentations\',...
                'caa_22\frontal\segmentations\',...
                'caa_22\occipital\segmentations\',...
                'caa_25\occipital\segmentations\',...
                'caa_26\occipital\segmentations\'};
    epvs_seg = {'caa17_occipital_THRESH-0.5_masked.mat',...
              'caa22-frontal_vessels-masked.mat',...
              'caa22-occipital_vessels-masked.mat',...
              'caa25-occipital_vessels-masked.mat',...
              'caa26-occipital_vessels-masked.mat'};
    % EPVS matrix crop coordinates
    crop_idcs = [1149,1549,1105,1705,1,596;...
                100,1512, 900,1300, 1,625;...
                200,1000, 200,1200, 1,353;...
                932,1432, 1084,1584,1,500;...
                700, 1429, 1000, 1581, 1, 592];
    % EPVS figure titles
    epvs_tit = {'CAA 17 Occipital','CAA 22 Frontal', 'CAA 22 Occipital',...
        'CAA 25 Occipital', 'CAA 26 Occipital'}; 
end

%% CAA22-Occipital: 
% Dynamically adjust the minimum branch length according to the length of
% the segmentation.

%%% File paths
caa22_occ_path = fullfile(dpath,'caa22/occipital/segmentations/');
caa22_occ_vasc_loops_fpath = 'caa22_occipital_vessels_crop_350_graph_data.mat';
caa22_occ_epvs_loops_fpath = 'caa22_occipital_epvs_crop_350_graph_data.mat';
caa22_occ_vasc_fpath = 'caa22_occipital_vessels_crop_350_rmloop_graph_data.mat';
caa22_occ_epvs_fpath = 'caa22_occipital_epvs_crop_350_rmloop_graph_data.mat';

%%% CAA 22 Occipital (with loops)
% Tissue Mask Crop
caa22_occ_mask_crop_350 = fullfile(caa22_occ_path, 'caa22_occipital_mask_edited_crop_350.tif');
caa22_occ_mask_crop_350 = TIFF2MAT(caa22_occ_mask_crop_350);
% EPVS Segmentation
caa22_occ_epvs_loops = load(fullfile(caa22_occ_path, caa22_occ_epvs_loops_fpath));
caa22_occ_epvs_angio = caa22_occ_epvs_loops.Data.angio;
caa22_occ_epvs_loops = caa22_occ_epvs_loops.Data;

%%% Adjust skeletonization minimum branch length
% Iterate over each segmentation
% Identify the length of the segmentation
% Adjust the skeletonization minimum branch length accordingly

% Identify groups of voxels
cc = bwconncomp(caa22_occ_epvs_angio);
% Measure the principal axes of the fitted ellipsoid of each group
s = regionprops3(cc, "PrincipalAxisLength");
% Separate into segmentation groups
pix_list = cc.PixelIdxList;
% Initialize matrix for storing skeletonization
sk = zeros(size(caa22_occ_epvs_angio));

% Iterate over each segmentation group
for ii = 1:length(pix_list)
    %%% Create segmentation matrix for group ii of voxels
    % Initialize a matrix for the segmentation group
    seg_mat = zeros(size(caa22_occ_epvs_angio));
    % Extract segmentation group pixel indices
    seg_idx = pix_list{ii};
    % Add ones to the segmentation matrix
    seg_mat(seg_idx) = 1;
    
    %%% initial approximation of minimum branch length 
    % Identify max principal length of ellipsoid of segmentation
    pmax = max(table2array(s(ii,:)));
    % Double the principal axis and subtract a small delta
    pmax = 2 .* pmax - 5;
    % Account for cases of segmentation one voxel in length
    pmax = max(1, floor(pmax - 5));
    
    %%% Generate skeleton for segment
    % Initialize matrix for storing the skeleton
    skel = [];
    % Continue decreasing branch length until skeleton is made
    while ~ismember(1,skel)
        % Adjust the skeletonization minimum branch length
        skel = bwskel(logical(seg_mat),'MinBranchLength',pmax);
        % Decrease pmax
        pmax = max(1, floor(pmax-5));
        if ~ismember(1,skel)
            fprintf('Empty on iteration %i\n',ii)
        end
    end
    % Add skeletonization to the skeleton matrix
    sk = sk + skel;
end


%%% Overlays: skeleton and segmentation (vessels and EPVS with loops)
volshow_overlay(sk, caa22_occ_epvs_angio,...
                'CAA22-Occi Vessel + Skeleton (adaptive branch length)')

%% CAA22-Frontal: Compare metrics before/after loop removal
% Use cropped regions of subjects CAA22-occipital and CAA22-frontal.
% Compare metrics (length density, branch density, tortuosity) from the
% skeleton of the EPVS and vascular segmentation

%%% File paths
caa22_front_path = fullfile(dpath,'caa22/frontal/segmentations/');
caa22_front_vasc_loops_fpath = 'caa22-frontal_vessels-masked_crop_350_graph_data.mat';
caa22_front_epvs_loops_fpath = 'caa22_frontal_epvs_min15vox_crop_350_graph_data.mat';
caa22_front_vasc_fpath = 'caa22-frontal_vessels-masked_crop_350_rmloop_graph_data.mat';
caa22_front_epvs_fpath = 'caa22_frontal_epvs_min15vox_crop_350_rmloop_graph_data.mat';

%%% CAA 22 Occipital (with loops)
% Tissue Mask Crop
caa22_front_mask_crop_350 = fullfile(caa22_front_path, 'caa22_frontal_mask_edited_crop_350.tif');
caa22_front_mask_crop_350 = TIFF2MAT(caa22_front_mask_crop_350);
% Vessel Segmentation
caa22_front_seg_loops = load(fullfile(caa22_front_path, caa22_front_vasc_loops_fpath));
caa22_front_seg_angio = caa22_front_seg_loops.Data.angio;
caa22_front_seg_loops = caa22_front_seg_loops.Data;
% EPVS Segmentation
caa22_front_epvs_loops = load(fullfile(caa22_front_path, caa22_front_epvs_loops_fpath));
caa22_front_epvs_angio = caa22_front_epvs_loops.Data.angio;
caa22_front_epvs_loops = caa22_front_epvs_loops.Data;
% Calculate metrics (vasculature)
caa22_front = struct();
caa22_front.ves.loops.bd = branch_density(caa22_front_seg_loops, caa22_front_mask_crop_350);
caa22_front.ves.loops.tort = calc_tortuosity(caa22_front_seg_loops);
caa22_front.ves.loops.mean_tort = mean(caa22_front.ves.loops.tort);
caa22_front.ves.loops.std_tort = std(caa22_front.ves.loops.tort);
caa22_front.ves.loops.ld = length_density(caa22_front_seg_loops, caa22_front_mask_crop_350);
% Calculate metrics (epvs)
caa22_front.epvs.loops.bd = branch_density(caa22_front_epvs_loops, caa22_front_mask_crop_350);
caa22_front.epvs.loops.tort = calc_tortuosity(caa22_front_epvs_loops);
caa22_front.epvs.loops.mean_tort = mean(caa22_front.epvs.loops.tort);
caa22_front.epvs.loops.std_tort = std(caa22_front.epvs.loops.tort);
caa22_front.epvs.loops.ld = length_density(caa22_front_epvs_loops, caa22_front_mask_crop_350);
% Overlays: skeleton and segmentation (vessels and EPVS with loops)
ves_skel = bwskel(caa22_front_seg_angio,'MinBranchLength',25);
volshow_overlay(ves_skel, caa22_front_seg_angio,'CAA22-Occi Vessel + Skeleton (min branch 25)')
epvs_skel = bwskel(caa22_front_epvs_angio,'MinBranchLength',25);
volshow_overlay(epvs_skel, caa22_front_epvs_angio,'CAA22-Occi Vessel + Skeleton (min branch 25)')
volshow_overlay(ves_skel(1:100,1:100,1:100),...
                caa22_front_seg_angio(1:100,1:100,1:100),...
                'CAA22-Occi Vessel + Skeleton (min branch 25)')

