%% Quality Assurance of CAA segmentations
% The purpose of this script is the following:
%   - CAA EPVS segmentation heuristics
%       - Overlay EPVS and segmentation in same ROI
%       - capture screenshot
%       - measure DICE score of EPVS
%           (manual EPVS segmentation vs. vascular segmentation from Eti)
%   - Generate Overlays of OCT, segmentation, & skeleton
%       - deep learning vascular segmentation (from Etienne)
%       - EPVS segmentation
%       - identify a large loop and a small loop
%   - Run graph/loop removal on EPVS segmentation
%   - Compare metrics before/after loop removal on a Sub-volume
%       - Take sub-volume
%       - generate graph before/after loop removal
%       - report metrics before/after
%       - branch density, length density, tortuosity
%   - Compare metrics between a simple & nested loop
%       - branch density, length density, tortuosity

clear; close all; clc;

%% Add top-level directory of code repository to path
% Start in current directory
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
% Set maximum number of threads equal to number of threads for script
NSLOTS = str2num(getenv('NSLOTS'));
maxNumCompThreads(NSLOTS);

%% Initialize data path for datasets
% Path to top-level directory
dpath = '/projectnb/npbssmic/ns/CAA/';
% Subfolder containing data
subdir = 'segmentations/';
% filename extension
ext = '.mat';

% Subject IDs
subids = {'caa6/frontal/', 'caa6/occipital/',...
         'caa17/occipital/',...
         'caa22/frontal/','caa22/occipital/',...
         'caa25/frontal/', 'caa25/occipital/',...
         'caa26/frontal/', 'caa26/occipital/'};

% Masked segmentations
seg_names = {'caa6-frontal_vessels-masked.mat',...
          'caa6-occipital_vessels-masked.mat',...
          'caa17_occipital_THRESH-0.5_masked.mat',...
          'caa22-frontal_vessels-masked.mat',...
          'caa22-occipital_vessels-masked.mat',...
          'caa25-frontal_vessels-masked.mat',...
          'caa25-occipital_vessels-masked.mat',...
          'caa26-frontal_vessels-masked.mat',...
          'caa26-occipital_vessels-masked.mat'};

% Skeleton corresponding to each segmentation/graph
skels = {'caa6-frontal_vessels-masked_skeleton.mat',...
          'caa6-occipital_vessels-masked_skeleton.mat',...
          'caa17_occipital_THRESH-0.5_masked_skeleton.mat',...
          'caa22-frontal_vessels-masked_skeleton.mat',...
          'caa25-frontal_vessels-masked_skeleton.mat',...
          'caa25-occipital_vessels-masked_skeleton.mat',...
          'caa26-frontal_vessels-masked_skeleton.mat',...
          'caa26-occipital_vessels-masked_skeleton.mat'};

% EPVS directory and filenames
epvs_paths = {'/caa17/occipital/segmentations/',...
              '/caa17/occipital/segmentations/',...
              '/caa22/frontal/segmentations/',...
              '/caa22/occipital/segmentations/',...
              '/caa25/occipital/segmentations/',...
              '/caa26/occipital/segmentations/'};

epvs_names = {'segmentation_07072023_crop.mgz',...
            'segmentation_07072023.mgz',...
            'EPVS_segmentation_03262024.mgz',...
            'EPVS_segmentation_02132024_registered.mgz',...
            'EPVS_mus_segmentation.mgz',...
            'EPVS_mus_segmentation.mgz'};

%% Load EPVS segmentations (.MGZ)
% Generate skeleton from segmentation
% Save both segmentation and skeleton as .MAT
%{
% Iterate over each sample
for ii = 1:length(epvs_names)
    % Define entire filepath
    fullpath = fullfile(dpath,epvs_paths{ii},epvs_names{ii});
    % Import the MGZ file
    epvs = MRIread(fullpath,0,1);
    epvs = logical(epvs.vol);
    % Convert segmentation to skeleton
    epvs_skel = bwskel(epvs);
    
    %%% Save the EPVS segmentation and skeleton
    % If the filename contains "crop" then append to output filename
    if contains(epvs_names{ii},'crop')
        epvs_out = fullfile(dpath,epvs_paths{ii},'epvs_crop.mat');
        skel_out = fullfile(dpath,epvs_paths{ii},'epvs_skel_crop.mat');
    else
        epvs_out = fullfile(dpath,epvs_paths{ii},'epvs.mat');
        skel_out = fullfile(dpath,epvs_paths{ii},'epvs_skel.mat');
    end
    save(epvs_out,'epvs','-v7.3');
    save(skel_out,'epvs_skel','-v7.3');
    clear epvs; clear epvs_skel;
end
%}

%% Overlay of EPVS crop + vasculature & DICE score


%% Load segmentations and skeletons
% Define entire filepath 
% fullpath = fullfile(dpath, subid, subdir);












