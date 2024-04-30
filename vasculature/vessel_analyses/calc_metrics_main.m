%% Main script for calling all the metric functions for n subjects. 
% Then, storing all of the metrics for each subject in a strict, within an array
% of structs with n elements (one for each subject)
% Functions being called: ave_length.m, branch_density.m,
% calc_tortuosity.m, length_density.m

%% Clear workspace & add top-level directory
clear; clc; close all;

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
    'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/vox_min_100'];

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

subid = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_7126',...
         'NC_6839',  'NC_6974',  'NC_8653',  'NC_21499', 'NC_301181'};

% Creating struct for storing all of the metrics for each subject
metrics = struct();

%% Iterate through subjects
for ii = 1:length(subid)
    % Specify the subject ID
    sub = subid{ii};
    fprintf('Calculating Metrics for Subject %s\n', sub);
    
    %%% Create struct for storing masks
    masks = struct();

    % Load all premade masks (tissue, gyri, sulci, GM, WM)
    mask_tiss = load(fullfile(dpath,sub,maskdir,'mask_tiss.mat'),'mask_tiss');
    mask_gyri = load(fullfile(dpath,sub,maskdir,'mask_gyri.mat'));
    mask_sulci = load(fullfile(dpath,sub,maskdir,'mask_sulci.mat'));
    mask_gm = load(fullfile(dpath,sub,maskdir,'mask_gm.mat'));
    mask_wm = load(fullfile(dpath,sub,maskdir,'mask_wm.mat'));
    mask_tiss = mask_tiss.mask_tiss;
    mask_gyri = mask_gyri.mask_gyri;
    mask_sulci = mask_sulci.mask_sulci;
    mask_gm = mask_gm.mask_gm;
    mask_wm = mask_wm.mask_wm;

    %%% Find the maximum depth to use for the masks
    zmin = min([size(mask_tiss,3), size(mask_gyri,2), size(mask_sulci,3),...
        size(mask_gm,3), size(mask_wm,3)]);

    %%% Truncate each volume to the maximum depths
    masks.tiss =     mask_tiss(:, :, 1:zmin);
    masks.wm =       mask_wm(:, :, 1:zmin);
    masks.gm =       mask_gm(:, :, 1:zmin);
    masks.sulci =    mask_sulci(:, :, 1:zmin);
    masks.gyri =     mask_gyri(:, :, 1:zmin);

    % Create remaining masks (WM & GM partitions of sulci/gyri)
    masks.gm_sulci = masks.gm .* masks.sulci;
    masks.wm_sulci = masks.wm .* masks.sulci;
    masks.gm_gyri = masks.gm .* masks.gyri;
    masks.wm_gyri = masks.wm .* masks.gyri;

    %% Iterate through each graph for each subject
    for j=1:9
        % Load the graph
        data = load(fullfile(dpath, sub, graphdir, graphs{j}));
        data = data.Data;
        
        % Load the respective mask for this iteration
        f = fields(masks);
        mask = masks.(f{j});
        
        % Call metric functions
        [um_len, ~] =   ave_length(data);
        lengths     =   data.Graph.segInfo.segLen_um;
        ld =            length_density(data, mask);
        branchden =     branch_density(data, mask);
        tort =          calc_tortuosity(data);
        fv =            fraction_vol(data, mask);
        
        % Store metrics for this subject and partition of the volume
        metrics.(sub).(f{j}).ave_length_um = um_len;
        metrics.(sub).(f{j}).lengths = lengths;
        metrics.(sub).(f{j}).length_density = ld;
        metrics.(sub).(f{j}).branch_density = branchden;
        metrics.(sub).(f{j}).tortuosity = tort;
        metrics.(sub).(f{j}).fraction_volume = fv;
    end

    fprintf('Finished Subject %s\n', sub);
end

%% Save output (maybe)

save(fullfile(mpath, 'metrics.mat'), 'metrics');