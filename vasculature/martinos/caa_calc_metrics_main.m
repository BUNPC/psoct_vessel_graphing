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
if ispc
    dpath = 'C:\Users\mhyman\boas_wang_lab\CAA\';
    % Metrics output path
    mpath = 'C:\Users\mhyman\boas_wang_lab\CAA\metrics\';
else
    dpath = '/projectnb/npbssmic/ns/CAA/';
    % Metrics output path
    mpath = '/projectnb/npbssmic/ns/CAA/metrics/';
end
   
% Directory containing masks
mdirs = {'caa6/frontal', 'caa6/occipital',...
        'caa17/occipital',...
        'caa22/frontal',...
        'caa25/frontal','caa25/occipital',...
        'caa26/frontal','caa26/occipital'};

% Name of subjects-region to save in struct
subs = strrep(mdirs,'/','_');

% Subdirectory containing the segmentation and graph
segdir = '/segmentations/';



% Graph structures to import and analyze
graphs = {'caa6-frontal_vessels-masked_crop_300_rmloop_graph_data.mat',...
          'caa6-occipital_vessels-masked_crop_300_rmloop_graph_data.mat',...
          'caa17_occipital_THRESH-0.5_crop_300_rmloop_graph_data.mat',...
          'caa22-frontal_vessels-masked_crop_300_rmloop_graph_data.mat',...
          'caa25-frontal_vessels-masked_crop_300_rmloop_graph_data.mat',...
          'caa25-occipital_vessels-masked_crop_300_rmloop_graph_data.mat',...
          'caa26-frontal_vessels-masked_crop_300_rmloop_graph_data.mat',...
          'caa26-occipital_vessels-masked_crop_300_rmloop_graph_data.mat'};

% Masks corresponding to each respective graph
masks = {'mask_crop_300.tif',...
        'mask_crop_300.tif',...
        'mask_crop_300.tif',...
        'mask_crop_300.tif',...
        'caa25_frontal_mask_edited_crop_300.tif',...
        'caa25_occipital_mask_edited_crop_300.tif',...
        'caa26_frontal_mask_edited_crop_300.tif',...
        'caa26_occipital_mask_crop_300.tif'};

% Creating struct for storing all of the metrics for each subject
metrics = struct();

% Change slashes if running on PC
if ispc
    mdirs = strrep(mdirs, '/','\');
    segdir = strrep(segdir, '/','\');
    mpath = strrep(mpath, '/','\');
end

%% Iterate through each segmentation graph
for ii=1:length(graphs)
    %%% Specify the subject and location
    sub = subs{ii};
    fprintf('Calculating Metrics for Subject %s\n', sub);
    
    %%% Create struct for storing masks
    % Load mask
    mask_path = fullfile(dpath, mdirs{ii}, segdir, masks{ii});
    [~,~,ext] = fileparts(masks{ii});
    if strcmp(ext,'.nii')
        mask = MRIread(mask_path, 0, 0);
        mask = logical(mask.vol);
    elseif strcmp(ext,'.tif')
        mask = TIFF2MAT(mask_path);
    end
    % Load the graph
    data = load(fullfile(dpath, mdirs{ii}, segdir, graphs{ii}));
    data = data.Data;
    
    % Call metric functions
    [um_len, ~] =   ave_length(data);
    lengths     =   data.Graph.segInfo.segLen_um;
    ld =            length_density(data, mask);
    branchden =     branch_density(data, mask);
    tort =          calc_tortuosity(data);
    fv =            fraction_vol(data, mask);

    % Calculate mean, median, mode of tortuosity
    tort_mean = mean(tort);
    tort_med = median(tort);
    tort_mode = mode(tort);
    
    % Store metrics for this subject and partition of the volume
    metrics.(sub).length_density = ld;
    metrics.(sub).branch_density = branchden;
    metrics.(sub).fraction_volume = fv;
    metrics.(sub).tortuosity_array = tort;
    metrics.(sub).tortuosity_mean = tort_mean;
    metrics.(sub).tortuosity_med = tort_med;
    metrics.(sub).tortuosity_mode = tort_mode;

    fprintf('Finished Subject %s\n', sub);
end

%% Save output
save(fullfile(mpath, 'caa_crop300_rmloops_metrics.mat'), 'metrics');