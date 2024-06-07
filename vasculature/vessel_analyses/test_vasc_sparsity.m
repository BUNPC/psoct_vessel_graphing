%% Test script for measuring vascular sparsity
% TODO:
%{
- 
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

%% Define directory paths

% Path to top-level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
   
% Volume and graph directories
voldir = '/dist_corrected/volume/ref/';
maskdir = 'dist_corrected/volume/ref/masks/';
graphdir = ['/dist_corrected/volume/combined_segs/' ...
            'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/vox_min_50/'];

% Metrics output path
mpath = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
    'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/vox_min_50/'];

% Graph structures to analyze
graphs = {'seg_refined_masked_rmloop_graph_data.mat',...
          'seg_gyri_refined_masked_rmloop_graph_data.mat',...
          'seg_sulci_refined_masked_rmloop_graph_data.mat',...
          'seg_gm_refined_masked_rmloop_graph_data.mat',...
          'seg_gm_gyri_refined_masked_rmloop_graph_data.mat',...
          'seg_gm_sulci_refined_masked_rmloop_graph_data.mat',...
          'seg_wm_refined_masked_rmloop_graph_data.mat',...
          'seg_wm_gyri_refined_masked_rmloop_graph_data.mat',...
          'seg_wm_sulci_refined_masked_rmloop_graph_data.mat'};

% Masks corresponding to each respective graph
masks = {'mask_tiss.mat','mask_gyri','mask_sulci','mask_gm','mask_wm'};

% IDs of each subject
subid = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_7126',...
         'NC_6839',  'NC_6974',  'NC_8653',  'NC_21499', 'NC_8095'};

%% Import the mask and segmentation
sub = subid{1};
mask_tiss = load(fullfile(dpath,sub,maskdir,'mask_tiss.mat'),'mask_tiss');
mask = mask_tiss.mask_tiss;
data = load(fullfile(dpath, sub, graphdir, graphs{1}));
data = data.Data;
angio = data.angio;

%% Calculate vascular sparsity
% Voxel dimensions (microns)
vox_dim = [12, 12, 15];
% Call function to calculate vascular sparsity
% vs = vasc_sparsity(angio, mask, vox_dim);


%%% Identify indices of non-vessel voxels in the masked OCT matrix
% Remove the vessel indices from tissue indices
tiss_idx = mask;
tiss_idx(angio) = [];
% Create array of non-vessel tissue indices (non-zero)
tiss_idx = find(tiss_idx);

%%% Convert vessel indices to subscripts
% Vessel indices (non-zero) in segmentation
vess_idx = find(angio);
% Convert vessel indices into subscript (y,x,z coordinate)
[vy, vx, vz] = ind2sub(size(angio),vess_idx);
vess = [vy, vx, vz];

%% Calculate difference between indices
% The euclidean distance function takes 0.6 seconds, and this will result
% in a computation time of 3000 days when operating on all non-vessel
% tissue voxels. The purpose of this section is to compute the minimum
% difference between indices, since these are all one-dimensional values.
% This code will identify the smallest difference in indices, and then it
% will compute the distance between the respective subcripts.

% Iterate over tissue indices
% tic
for ii = 1:size(tiss_idx,2)
    tic
    t = tiss_idx(ii);
    % Find difference between tissue index and all vessel indices
    d = zeros(size(vess_idx,1),1);
    for j=1:length(vess_idx)
        d(j) = abs(t-vess_idx(j));
    end
    % Find the smallest difference
    dmin = min(d);
    toc
end
% toc

%% Find distance between non-vessel tissue voxel and all vessel voxels
% Initialize array to store the minimum distance between tiss/vessel
dmin = zeros(1,size(tiss_idx,2));
% Convert tissue indices to subscripts
[ty, tx, tz] = ind2sub(size(angio),tiss_idx);
tiss_sub = [ty',tx',tz'];

for ii = 1:size(tiss_idx,2)
    % Find euclidean distance between tissue and vessel
    t = tiss_sub(ii,:);
    d = zeros(size(vess,1),1);
    tic
    for j = 1:length(vess)
        % Take the jth element of the vessels
        v = vess(j,:);
        % Compute distance between vessel and tissue
        d(j) = sqrt((v(1)-t(1)).^2 + (v(2)-t(2)).^2 + (v(3)-t(3)).^2);        
    end
    toc
    % Find the minimum distance for this tissue index
    dmin(ii) = min(d);
end




