%% Main script for measuring vascular sparsity
% Iterate over each subject and tissue region for each subject.

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
            'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];

% Metrics output path
mpath = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
    'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];

% Graph structures to analyze
graphs = {'seg_refined_masked_rmloop_graph_data.mat',...
          'seg_gyri_refined_masked_rmloop_graph_data.mat',...
          'seg_sulci_refined_masked_rmloop_graph_data.mat',...
          'seg_gm_refined_masked_rmloop_graph_data.mat',...
          'seg_wm_refined_masked_rmloop_graph_data.mat',...
          'seg_gm_sulci_refined_masked_rmloop_graph_data.mat',...
          'seg_wm_sulci_refined_masked_rmloop_graph_data.mat',...
          'seg_gm_gyri_refined_masked_rmloop_graph_data.mat',...
          'seg_wm_gyri_refined_masked_rmloop_graph_data.mat'};

% Tissue Regions
regions = {'tiss','gyri','sulci','gm','wm','gm_sulci','wm_sulci',...
            'gm_gyri','wm_gyri'};

% Masks corresponding to each respective graph
masks = {'mask_tiss.mat','mask_gyri.mat','mask_sulci.mat',...
        'mask_gm.mat','mask_wm.mat'};

% IDs of each subject
subid = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_7126',...
         'NC_6839',  'NC_6974',  'NC_8653',  'NC_21499', 'NC_8095'};

% Voxel dimensions (microns)
vox_dim = [12, 12, 15];

%% Iterate over subjects and regions - calculate vascular sparsity

% Struct to store the vascular sparsity values
vs = struct();

% Iterate over subjects
for sidx = 1 : length(subid)
    % Retrieve subject ID
    sub = subid{sidx};
    
    %%% Generate masks
    % There are only masks for the entire tissue volume, gyri, sulci, GM,
    % and WM. This section generates teh masks gm_sulci, wm_sulci, gm_gyri,
    % wm_gyri.
    % 
    % Load the five masks that have already been generated
    mstruct = struct();
    for m = 1:5
        mask = load(fullfile(dpath,sub,maskdir,masks{m}));
        f = fields(mask);
        mstruct.(regions{m}) = mask.(f{1});
    end
    % Generate the remaining masks
    mstruct.(regions{6}) = logical(mstruct.gm .* mstruct.sulci);
    mstruct.(regions{7}) = logical(mstruct.wm .* mstruct.sulci);
    mstruct.(regions{8}) = logical(mstruct.gm .* mstruct.gyri);
    mstruct.(regions{9}) = logical(mstruct.wm .* mstruct.gyri);
    
    % Iterate over sub-regions
    for r = 1:length(graphs)
        % Import the respective mask
        mask = mstruct.(regions{r});

        % Import the vessel segmentation (angio)
        data = load(fullfile(dpath, sub, graphdir, graphs{sidx}));
        data = data.Data;
        angio = data.angio;

        % Calculate vascular sparsity for subject and sub-region
        vs.(sub).(regions{r}) = vasc_sparsity(angio, mask, vox_dim);

        % Print to console
        fprintf('Completed region "%s" of subject %s\n',regions{r}, sub)
    end
    fprintf('FINISHED SUBJECT %s\n\n',sub)
end

% Save the vascular sparsity output
fout = fullfile(mpath, 'metrics_vasc_sparsity.mat');
save(fout,'vs','-v7.3');