%% Vessel metrics - statistical analyses
%{
This script is for analyzing the length density in the dataset:
/projectnb/npbssmic/ns/Ann_Mckee_samples_10T
%}
clear; clc; close all;

%% Add top-level directory of code repository to path
% This allows Matlab to find the functions in the project folders

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

%% Data Directories and filenames
% Path to top-level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';

% Additional subfolder when applying minimum voxel count
% vox_min = 'vox_min_50';
vox_min = false;

if vox_min
    % Volume and graph directories
    voldir = '/dist_corrected/volume/ref/';
    graphdir = ['/dist_corrected/volume/combined_segs/' ...
                'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];
    graphdir = fullfile(graphdir, vox_min);
    % Metrics output path
    mpath = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
        'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];
    mpath = fullfile(mpath, vox_min);
else
    % Volume and graph directories
    voldir = '/dist_corrected/volume/ref/';
    graphdir = ['/dist_corrected/volume/combined_segs/' ...
                'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];
    
    % Metrics output path
    mpath = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
        'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];
end

%% Import the metrics struct from the last script
% Load the metrics struct and extract subject IDs
metrics = load(fullfile(mpath, 'metrics.mat'));
metrics = metrics.metrics;
subids = fields(metrics);
metrics_out = fullfile(mpath, 'metrics.mat');

% Output filename for saving the p-value table
ptable_out = 'p_value_table.xls';

% Voxel dimensions (microns) and volume (cubic micron)
vox_dim = [12, 12, 15];
vox_vol = vox_dim(1) .* vox_dim(2) .* vox_dim(3);

% Region of brain (excluding ratios)
regions = {'tiss','gyri','sulci','gm','wm','gm_sulci','wm_sulci',...
            'gm_gyri','wm_gyri'};
% Regions of brain (including ratios)
all_regions = {'tiss','gyri','sulci','gm','wm','gm_sulci','wm_sulci',...
               'gm_gyri','wm_gyri','sulci_gyri','wm_sulci_gyri',...
               'gm_sulci_gyri'};

%% Calculate N tortuosity outliers for each subject
% Calculate the outlier threshold from the healthy control group
% Then, calculate the number of outliers for each group (for each tissue
% region)

%%% Calculate outlier threshold based on the healthy control group
for ii = 1:length(regions)
    % Combine tortuosity measurements from all normal controls
    [~, ~, nc] = organize_metrics(metrics, subids,regions{ii}, 'tortuosity');
    
    % Calculate outlier cutoff = 3*MAD (scaled). This is taken from the
    % Mathworks website for calculating the threshold
    c = -1 / (sqrt(2)*erfcinv(3/2));
    mad3 = 3*c*mad(nc,1);
    cutoff = median(nc) + mad3;

    %%% Calculate number of outliers for each subject
    for j=1:length(subids)
        % Load tortuosity for each subject/region
        t = metrics.(subids{j}).(regions{ii}).tortuosity;
        n = sum(t > cutoff);
        % Save the number of outliers to struct
        metrics.(subids{j}).(regions{ii}).tort_outliers = n;
    end
end

%% Reorganize data and generate barcharts
% This script will reorganize the metrics struct by taking the region from
% each subject and combining them into a single substruct (e.g. taking the
% "tiss" region from each subject and creating a substruct named
% "metrics.tiss". This section can be skipped if it was already run.

% Metric
params = {'length_density','branch_density','fraction_volume',...
    'tort_outliers','tortuosity'};
% Title of each bar chart
titles = {'Length Density','Branch Density','Volume Fraction',...
    'N Tortuosity Outliers'};
% Region titles
r_titles = {'Tissue', 'Gyri', 'Sulci', 'GM', 'WM',...
            'GM Sulci', 'WM Sulci', 'GM Gyri', 'WM Gyri'};
% Y-axis labels
ylabels = {'Length Density (\mum^-^2)','Branch Density (\mum^-^3)',...
            'Volume Fraction (a.u.)','Count'};
% Name of file to save
plot_names = {'length_density','branch_density','volume_fraction'...
    'tortuosity_outliers'};

% Iterate over each tissue region
for ii = 1:length(regions)
    % Iterate over metrics
    for j = 1:length(params)
        % Retrieve the metrics for this region/parameter
        [ad, cte, nc] = organize_metrics(metrics, subids,...
                                         regions{ii}, params{j});
        % Save parameter by the group to fascilitate statistical analyses
        % in a later step.
        metrics.(regions{ii}).(params{j}).ad = ad;
        metrics.(regions{ii}).(params{j}).cte = cte;
        metrics.(regions{ii}).(params{j}).nc = nc;

        % Create and save the bar chart (except for tortuosity
        if j~=5
            % Filename to save
            plot_name = strcat(regions{ii},'_',plot_names{j});
            % Title for bar chart
            tstr = append(r_titles{ii},' - ', titles{j});
            % Create bar chart
            bar_chart(ad, cte, nc, subids, mpath, tstr,'',...
                      ylabels{j}, plot_name);
            close;
        end
    end
end

%% Calculate ratio of sulci/gyri
% This cannot be calculated for the tortuosity since it's an array of
% values.

% Cell array of parameters to compute ratio 
ratio_params = {'length_density','branch_density','fraction_volume',...
    'tort_outliers'};
groups = {'ad','cte','nc'};

% Iterate over each parameter
for ii = 1:length(ratio_params)
    % Iterate over groups
    for j = 1:length(groups)
        %%% Retrieve the parameter for each region
        sulci = metrics.sulci.(ratio_params{ii}).(groups{j});
        gyri = metrics.gyri.(ratio_params{ii}).(groups{j});
        wm_sulci = metrics.wm_sulci.(ratio_params{ii}).(groups{j});
        gm_sulci = metrics.gm_sulci.(ratio_params{ii}).(groups{j});
        wm_gyri = metrics.wm_gyri.(ratio_params{ii}).(groups{j});
        gm_gyri = metrics.gm_gyri.(ratio_params{ii}).(groups{j});
        
        %%% Compute the ratios
        % sulci / gyri ratios
        sulci_gyri = sulci ./ gyri;
        % WM sulci / gyri ratios
        wm_sulci_gyri = wm_sulci ./ wm_gyri;
        % GM sulci / gyri ratios
        gm_sulci_gyri = gm_sulci ./ gm_gyri;

        %%% Save to the struct
        metrics.sulci_gyri.(ratio_params{ii}).(groups{j}) = sulci_gyri;
        metrics.wm_sulci_gyri.(ratio_params{ii}).(groups{j}) = wm_sulci_gyri;
        metrics.gm_sulci_gyri.(ratio_params{ii}).(groups{j}) = gm_sulci_gyri;
    end
end

%% Statistical Hypothesis Testing
% Trend threshold
trend = 0.10;
% Significant Difference threshold
alpha = 0.05;
% Calculate stats
pstats = metrics_stats(metrics, all_regions, params, alpha,...
                        trend, mpath, ptable_out);
% Save statistics
stats_fout = fullfile(mpath, 'stats.mat');
save(stats_fout, 'pstats','-v7.3');
%% Generate violin plots of tortuosity
% This section will iterate over each vascular metric and generate a 1xm
% cell array, where each column contains the vascular metric values for a
% particular brain region. Then, this section will generate a violin plot
% for each vascular metric, where the x-axis will represent different brain
% regions. For each brain region, the AD, CTE, NC will be grouped together

% Y-axis labels
ylabels = {'Length Density (\mum^-^2)','Branch Density (\mum^-^3)',...
            'Volume Fraction (a.u.)','Tortuosity (a.u.)'};
% X-axis labels for groups
xlabels = {'Volume','Gyri','Sulci','GM','WM','GM Sulci','WM Sulci',...
               'GM Gyri','WM Gyri'};
% Regions of the brain for calculating ratios
ratio_regions = {'sulci_gyri','wm_sulci_gyri','gm_sulci_gyri'};
ratio_ylabels = {'Length Density (\mum^-^2)','Branch Density (\mum^-^3)',...
            'Volume Fraction (a.u.)'};
ratio_xlabels = {'Sulci / Gyri', 'WM Sulci / Gyri', 'GM Sulci / Gyri'};

%%% Iterate over vascular metrics
for ii = 1:length(params)   
    % Initialize array to store the values
    vmetric = [];

    % Store disease state index (AD, CTE, HC)
    comp_idx = {};

    % Store brain region index ()
    region_idx = {};

    %% Iterate over brain regions (exclude ratios)
    for j = 1:length(regions)
        % Vertically concatenate the [AD; CTE; NC] into single vert array
        ad = metrics.(regions{j}).(params{ii}).ad;
        cte = metrics.(regions{j}).(params{ii}).cte;
        nc = metrics.(regions{j}).(params{ii}).nc;
        vmetric = vertcat(vmetric, ad, cte, nc);
        n_samples = length(ad) + length(cte) + length(nc);

        % Initialize labels for comparisons within each group (AD, CTE, HC)
        ad = repmat({'AD'},[length(ad),1]);
        cte = repmat({'CTE'},[length(cte),1]);
        nc = repmat({'HC'},[length(nc),1]);
        comp_idx = vertcat(comp_idx, ad, cte, nc);

        % Initialize label for comparisons between brain regions
        region_idx = vertcat(region_idx,...
            repmat(cellstr(xlabels{j}),[n_samples,1]));
    end

    %%% Initialize table for violin plots
    vtable = table(region_idx, comp_idx, vmetric);

    %%% Call function to generate violin plot
    % Group each AD/HC/CTE comparison by tissue
    grpandplot(vtable,"vmetric", yTitle = ylabels{ii},...
        xFactor="region_idx", cFactor = "comp_idx", xOrder = xlabels,...
        showXLine = true, showVln = true, showBox = false, showMean=true,...
        showPnt = false, showNum = false, numYPos = 500, pntSize=5);
    set(gca, 'FontSize', 30)
    % Save output
    fout = fullfile(mpath,append(params{ii},'_violin'));
    saveas(gcf,fout,'png');
    close;
    %% Iterate over brain regions (only ratios)
    if strcmp(params{ii}, 'tortuosity')
        continue
    else
        % Initialize array to store the values
        vmetric = [];
        % Store disease state index (AD, CTE, HC)
        comp_idx = {};
        % Store brain region index ()
        region_idx = {};
        
        % Iterate over ratio regions
        for j = 1:length(ratio_regions)
            % Vertically concatenate [AD; CTE; NC] into single vert array
            ad = metrics.(ratio_regions{j}).(params{ii}).ad;
            cte = metrics.(ratio_regions{j}).(params{ii}).cte;
            nc = metrics.(ratio_regions{j}).(params{ii}).nc;
            vmetric = vertcat(vmetric, ad, cte, nc);
            n_samples = length(ad) + length(cte) + length(nc);
    
            % Initialize labels for comparisons within group (AD, CTE, HC)
            ad = repmat({'AD'},[length(ad),1]);
            cte = repmat({'CTE'},[length(cte),1]);
            nc = repmat({'HC'},[length(nc),1]);
            comp_idx = vertcat(comp_idx, ad, cte, nc);
    
            % Initialize label for comparisons between brain regions
            region_idx = vertcat(region_idx,...
                repmat(cellstr(ratio_xlabels{j}),[n_samples,1]));
        end
    
        %%% Initialize table for violin plots
        vtable = table(region_idx, comp_idx, vmetric);
    
        %%% Call function to generate violin plot
        % Group each AD/HC/CTE comparison by tissue
        grpandplot(vtable,"vmetric", yTitle=ratio_ylabels{ii},...
            xFactor="region_idx", cFactor="comp_idx",...
            xOrder = ratio_xlabels, showXLine = true, showVln = true,...
            showBox = false, showMean=true,...
            showPnt = false, showNum = false, numYPos = 500, pntSize=5);
        set(gca, 'FontSize', 30)
        % Save output
        fout = fullfile(mpath,append(params{ii},'_ratio_violin'));
        saveas(gcf,fout,'png');
        close;
    end

end