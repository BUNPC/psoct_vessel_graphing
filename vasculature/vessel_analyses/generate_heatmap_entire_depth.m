%% Divide volume into grid and calculate metrics for each cube
% Generate heatmaps centered about specific depths in the OCT volumes.
% The purpose of this is to correlate the vascular metrics with the
% pathology staining of deposits of amyloid beta and phosphorylated tau.
% The output will be "heatmap_ab_ptau_[ROI size in microns].mat"

%% Add top-level directory of code repository to path
clear; clc; close all;
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
% Set maximum number of threads equal to number of threads for script
NSLOTS = str2num(getenv('NSLOTS'));
maxNumCompThreads(NSLOTS);

%% Initialize directories, filenames, parameters

%%% Directories 
% Metrics output path
mpath = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
    'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];

%%% Subvolume parameters
% Isotropic cube length (microns)
cube_side = 1000;
% Boolean for Minimum intensity projection of tissue mask
% true = minimum intensity projection, false = maximum intensity projection
min_ip = false;
% Size of each voxel (microns)
vox = [12, 12, 15];
% Whether to plot non-normalized heatmaps for each depth
viz_individual = false;

%%% Vascular metrics (volume fraction, length density, branch density,
% tortuosity, diameter)
metrics = {'vf','ld','bd','tort','diam'};

%% Vascular heatmap across entire heatmap depth
% Create a heatmap across the entire z-axis. Fig. 2 in manuscript. Exclude
% ROIs that are devoid of vasculature

% Load heatmap
hm = append('heatmap_',num2str(cube_side),'.mat');
hm = fullfile(mpath,hm);
hm = load(hm);
hm = hm.heatmap;
% Extract subject IDs from heatmap
subid = fields(hm);
% Initialize struct to store heatmaps across depth
hm_stack = struct();
% Initialize the maximum values for x,y dimensions. These are used for
% scaling all of the heatmap figures such that they are all on the same
% scale.
ymax = 1;
xmax = 1;

% Take MIP of each subject
for ii = 1:length(subid)
    % Take MIP of tissue mask
    sub = subid{ii};
    hm_stack.(sub).mask = max(hm.(sub).mask,[],3);
    % Identify the maximum dimensions of x,y
    ymax = max([ymax,size(hm_stack.(sub).mask,1)]);
    xmax = max([xmax,size(hm_stack.(sub).mask,2)]);
    
    %%% Take average of heatmap across z-axis   
    for j = 1:length(metrics)
        % Set the zero elements to NaN so that they're excluded from mean
        tmp = hm.(sub).(metrics{j});
        tmp(tmp==0) = NaN;
        % Take mean across z-dimension, exclude NaN
        hm_stack.(sub).(metrics{j}) = mean(tmp,3,'omitnan');
    end

    fprintf('FINISHED HEATMAP AVERAGE FOR SUBJECT %s\n',sub)
end

% Save the heatmap struct
heat_out = append('heatmap_entire_depth_',num2str(cube_side),'.mat');
heat_out = fullfile(mpath, heat_out);
save(heat_out,'hm_stack','-v7.3');

%% Calculate limits for normalized heatmaps
% Iterate over each metric, choose one depth from each subject, normalize
% the colorbar across all subjects for this metric.

% Load Heat map
heat_out = append('heatmap_entire_depth_',num2str(cube_side),'.mat');
heat_out = fullfile(mpath, heat_out);
heatmap = load(heat_out);
heatmap = heatmap.hm_stack;
subid = fields(heatmap);

%%% Identify minimum / maximum of each metric
% Set the maximum to be the 95th percentile of each metric. This will
% scale the colorbar to account for outliers.
vf = [];
ld = [];
bd = [];
tr = [];
dm = [];
for ii = 1:length(subid)
    vf = [vf, reshape(heatmap.(subid{ii}).vf(:),1,[])];
    ld = [ld, reshape(heatmap.(subid{ii}).ld(:),1,[])];
    bd = [bd, reshape(heatmap.(subid{ii}).bd(:),1,[])];
    tr = [tr, reshape(heatmap.(subid{ii}).tort(:),1,[])];
    dm = [dm, reshape(heatmap.(subid{ii}).diam(:),1,[])];
end
% Minimum of each metric
vf_min = min(vf);
ld_min = min(ld);
bd_min = min(bd);
tr_min = 1;
dm_min = min(dm);
% Maximum = 90th percentile of each metric
vf_max = prctile(vf,95);
ld_max = prctile(ld,95);
bd_max = prctile(bd,95);
tr_max = prctile(tr,95);
dm_max = prctile(dm,95);

%% Generate normalized heatmaps
% Variable for whether or not to invert the heatmap
flip_cbar = 0; 
% Iterate over each subject
for ii = 1:length(subid)
    %%% Load the heatmap for the subject
    sub = subid{ii};
    heatmap_vf = heatmap.(sub).vf;
    heatmap_ld = heatmap.(sub).ld;
    heatmap_bd = heatmap.(sub).bd;
    heatmap_tr = heatmap.(sub).tort;
    heatmap_dm = heatmap.(sub).diam;
    masks = heatmap.(sub).mask;
    
    %%% Rotate and transpose subjects so that they all align
    % Rotate subject AD_20969
    if strcmp(sub,'AD_20969')
        heatmap_vf = flip(heatmap_vf',2);
        heatmap_ld = flip(heatmap_ld',2);
        heatmap_bd = flip(heatmap_bd',2);
        heatmap_tr = flip(heatmap_tr',2);
        heatmap_dm = flip(heatmap_dm',2);
        masks = flip(masks',2);
    end
    if strcmp(sub,'CTE_6912')||strcmp(sub,'NC_6974')||strcmp(sub,'NC_21499')
        heatmap_vf = flip(heatmap_vf,1);
        heatmap_ld = flip(heatmap_ld,1);
        heatmap_bd = flip(heatmap_bd,1);
        heatmap_tr = flip(heatmap_tr,1);
        heatmap_dm = flip(heatmap_dm,1);
        masks = flip(masks,1);
    end

    %%% Output filepath for figures
    roi_dir = strcat('entire_depth_',num2str(cube_side));
    heatmap_dir = fullfile(mpath,'heatmaps',sub,roi_dir);
    if ~isfolder(heatmap_dir)
        mkdir(heatmap_dir);
    end

    %%% Plot heatmaps at each depth
    plot_save_heatmap([], heatmap_vf, flip_cbar, [vf_min, vf_max],...
        [ymax,xmax],masks,'Volume Fraction','(a.u.)',...
        heatmap_dir,'rescaled_heatmap_vf')
    plot_save_heatmap([], heatmap_ld, flip_cbar, [ld_min, ld_max],...
        [ymax,xmax],masks,'Length Density','Length (\mu) / Volume (\mu^3)',...
        heatmap_dir,'rescaled_heatmap_ld')
    plot_save_heatmap([], heatmap_bd, flip_cbar, [bd_min, bd_max],...
        [ymax,xmax],masks,'Branch Density','Branches / Volume (\mu^3)',...
        heatmap_dir,'rescaled_heatmap_bd')
    plot_save_heatmap([], heatmap_tr, flip_cbar, [tr_min, tr_max],...
        [ymax,xmax],masks,'Tortuosity','(a.u.)',...
        heatmap_dir,'rescaled_heatmap_tr')
    plot_save_heatmap([], heatmap_dm, flip_cbar, [dm_min, dm_max],...
        [ymax,xmax],masks,'Diameter','(\mum)',...
        heatmap_dir,'rescaled_heatmap_dm')
end


%% Plot and save the heat maps
function plot_save_heatmap(Ndepths, heatmaps, flip_cbar, colorbar_range,...
    max_dim, masks, tstr, cbar_label, dpath, fname)
% PLOT_SAVE_HEATMAP use imagesc and set background = 0
% INPUT
%   Ndepths (int): number of depths in z dimension
%   heatmaps (double matrix): heatmaps of vascular metric
%   flip_cbar (logical): reverse the direction of the colorbar
%   colorbar_range (double array): [min, max]
%   max_dim (double array): [ymax, xmax] the maximum dimensions of all
%                           heatmap matrices for the y- and x-axes.
%   masks (double): tissue mask (1=tissue, 0=other)
%   tstr (string): figure title
%   cbar_label (string): colorbar label
%   dpath (string): data directory path
%   fname (string): name of figure to save

%%% Set the number of depths to iterate for each heatmap
% If the number of depths is not specified, then set it equal to the number
% of z dimensions.
if isempty(Ndepths)
    Ndepths = size(heatmaps,3);
end
% Set fontsize for the heatmap figure
fontsize = 20;
% Set the maximum dimensions for x-y axes
ymax = max_dim(1);
xmax = max_dim(2);

%%% Iterate over frames in z dimension
for d = 1:Ndepths
    %%% Heatmap of j_th frame from the length density
    % If there are multiple heatmaps in the matrix
    if size(heatmaps,3) > 1
        heatmap = heatmaps(:,:,d);
    % Here it is just a single frame of a heatmap
    else
        heatmap = heatmaps;
    end
    % Initialize heatmap
    h = imagesc(heatmap);
    
    %%% Scale the x- and y-axes according to the max dimensions
    % Calculate the ratio of the current axes to the max axes
    yratio = size(heatmap,1) ./ ymax;
    xratio = size(heatmap,2) ./ xmax;
    set(gcf,'Units','Normalized','OuterPosition',[0,0,xratio,yratio])

    %%% Initialize colorbar
    % If the colorbar_range is passed in, then extract min & max
    if ~isempty(colorbar_range)
        cmap_min = colorbar_range(1);
        cmap_max = colorbar_range(2);
    % Otherwise, set limits from the current heatmap
    else
        % Min = Lowest value also greater than zero
        cmap_min = min(heatmap(heatmap(:)>0));
        % Max = Find the 95th percentile for upper limit
        cmap_max = prctile(heatmap(:), 95);
    end
    % Initialize the colormap limits
    cmap = jet(256);
    clim(gca, [cmap_min, cmap_max]);
    % Initialize colormap and colorbar
    if flip_cbar
        colormap(flipud(cmap));
    else
        colormap(cmap);
    end
    c = colorbar;

    %%% Apply tissue mask to the heatmap to remove background
    alpha_mask = double(masks(:,:,d));
    set(h, 'AlphaData', alpha_mask);

    %%% Configure figure parameters
    % Update title string with specific pathology
    pathology = {'A-Beta','p-tau'};
    if size(heatmaps,3) > 1
        title_str = append(tstr, ' ', pathology{d});
    else
        title_str = tstr;
    end
    title(title_str,'Interpreter','none');
    set(gca, 'FontSize', fontsize);
    % Label the colorbar    
    c.Label.String = cbar_label;
    % Offset colorbar label to the right of colorbar
    c.Label.Position = [10 (cmap_max - (cmap_max-cmap_min)/2)];
    c.Label.Rotation = 270;
    % Increase fontsize of colorbar
    c.FontSize = 20;
    
    %%% Save figure as PNG
    % If there are multiple heatmaps in the matrix, save vascular heatmap
    if size(heatmaps,3) > 1
        fout = append(fname, '_', pathology{d});
    % Otherwise, save the pathology heatmap
    else
        fout = fname;
    end
    % If the colorbar is reversed then add suffix to filename
    if flip_cbar
        fout = append(fout, '_flip_cbar');
    end
    % Save figure as PNG
    fout = fullfile(dpath, fout);
    saveas(gca, fout,'png');
    pause(0.1)
    close;
end
end