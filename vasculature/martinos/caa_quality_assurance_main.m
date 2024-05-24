%% Main file for calling segmentation functions
% Author: Mack Hyman
% Date Created: April 24, 2024
%
% Detailed Description
%{
1) 3D overlay of EPVS segmentation
    - EPVS segmentation + epvs_skeleton
    - EPVS segmentation + vascular segmentation
    - DICE score (EPVS vs. vascular)
2) 2D overlays
    - OCT, vessel segmentation, vessel epvs_skeleton
    - OCT, EPVS segmentation, EPVS epvs_skeleton
3) Compare metrics before/after loop removal
    - length density, branch density, tortuosity
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

%% Create 3D overlays (EPVS + epvs_skel, EPVS + vasculature)
%{
for ii = 3:length(epvs_dir)
    %%% EPVS + epvs_skeleton
    % Load the EPVS
    epvs = fullfile(dpath, epvs_dir{ii},'epvs.mat');
    epvs = load(epvs);
    epvs = epvs.epvs;
    % Crop the EPVS
    y0 = crop_idcs(ii,1);
    yf = crop_idcs(ii,2);
    x0 = crop_idcs(ii,3);
    xf = crop_idcs(ii,4);
    z0 = crop_idcs(ii,5);
    zf = crop_idcs(ii,6);
    epvs_crop = epvs(y0:yf, x0:xf, z0:zf);

    % % Load the epvs_skeleton
    % epvs_skel = fullfile(dpath, epvs_dir{ii},'epvs_skel.mat');
    % epvs_skel = load(epvs_skel);
    % epvs_skel = epvs_skel.epvs_skel;
    % % crop the epvs_skeleton
    % epvs_skel = epvs_skel(y0:yf, x0:xf, z0:zf);
    % 
    % % 3D overlay
    % volshow_overlay(epvs_skel, epvs, epvs_tit{ii});

    %%% EPVS + vasculature   
    % Load vasculature
    seg = fullfile(dpath, epvs_dir{ii},epvs_seg{ii});
    seg = load(seg);
    seg = seg.seg;
    if size(seg,1) == size(epvs,2)
        seg = permute(seg,[2,1,3]);
    end
    % crop vasculature
    seg_crop = seg(y0:yf, x0:xf, z0:zf);
    % mask vasculature with the EPVS
    seg_crop(~epvs_crop) = 0;
    % 3D overlay
    % volshow_overlay(seg_crop, epvs_crop, epvs_tit{ii});

    %%% DICE score
    if ii == 1
        epvs = epvs(:,:,1:end-4);
    end
    seg(~epvs) = 0;
    seg = logical(seg);
    d = dice(seg, epvs);
    sprintf('DICE score for %s = %f',epvs_tit{ii}, d)
end
%}

%% 2D overlays

%%% OCT, vessel segmentation, vessel skeleton
for ii=3:length(subids)
    % import segmentation
    seg = load(fullfile(dpath, subids{ii}, subdir, seg_names{ii}));
    seg = logical(seg.seg);

    % skeletonize segmentation
    % seg_skel = bwskel(seg);

    % import OCT
    vol = fullfile(dpath, subids{ii}, subdir, mus_name);
    if ~isfile(vol)
        vol = strcat(vol,'.gz');
    end
    vol = MRIread(vol,0,1); vol = vol.vol;
    vol = rescale(vol,0,1);
    vol = im2uint8(vol);
    
    % overlay OCT + segmentation
    if size(vol,1) == size(seg,2)
        vol = permute(vol,[2,1,3]);
    end
    % Truncate the volume for CAA 17 occipital
    if strcmp('caa_17/occipital/',subids{ii})
        vol = vol(:,:,1:end-4);
    end
    fout = fullfile(dpath, subids{ii}, subdir, 'oct_seg_overlay.tif');
    overlay_vol_seg(vol, seg, 'magenta', fout, false)

end

%%% OCT, EPVS segmentation, EPVS epvs_skeleton
for ii = 1:length(epvs_dir)
    %%% EPVS + epvs_skeleton
    % Load the EPVS
    epvs = fullfile(dpath, epvs_dir{ii},'epvs.mat');
    epvs = load(epvs);
    epvs = epvs.epvs;

    % skeletonize EPVS segmentation
    epvs_skel = bwskel(seg);

    % import OCT
    % TODO: update directory path to match the indexing
    vol = load(fullfile(dpath, subids{ii}, subdir, vol_names{ii}));
    
    %%% overlays
    % overlay OCT + segmentation
    fout = fullfile(dpath, subids{ii}, subdir, 'oct_epvs.tif');
    overlay_vol_seg(vol, epvs, 'magenta', fout, false)
    % overlay OCT + skeleton
    fout = fullfile(dpath, subids{ii}, subdir, 'oct_epvs_skel.tif');
    overlay_vol_seg(vol, epvs_skel, 'magenta', fout, false)
end


%% Compare metrics before/after loop removal
% length density, branch density, tortuosity

%%% metrics of segmentation with loops
% Load the segmentation of the cropped segmentation
seg = fullfile(dpath, subids{1}, subdir, 'caa6-frontal_vessels-masked_crop_300.tif');
seg = TIFF2MAT(seg);
seg = logical(seg);

% Create graph and skeleton (with loops)
fout = fullfile(dpath, subids{1}, subdir);
fname_out = 'caa6_frontal_crop';
viz = false;
rm_loop = false;
seg_graph_init(seg,[20,20,20],fout,fname_out,viz,rm_loop);

% Calculate metrics

% Create graph and skeleton (without loops)
[graph, skel] = seg_to_graph(seg, [20,20,20]);
[nodes, edges] = rm_loops(graph.nodes, graph.edges, seg, 6, 0.99, 1, false);
graph.nodes = nodes;
graph.edges = edges;
Data = init_graph(graph);













