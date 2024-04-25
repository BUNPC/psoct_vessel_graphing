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

%% Data Directories
% Path to top-level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
   
% Volume and graph directories
voldir = '/dist_corrected/volume/ref/';
graphdir = ['/dist_corrected/volume/combined_segs/' ...
            'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];

% Metrics output path
mpath = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
    'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];

% Load the metrics struct and extract subject IDs
metrics = load(fullfile(mpath, 'metrics.mat'));
metrics = metrics.metrics;
subids = fields(metrics);


% Voxel dimensions (microns) and volume (cubic micron)
vox_dim = [12, 12, 15];
vox_vol = vox_dim(1) .* vox_dim(2) .* vox_dim(3);

%% Reorganize data and generate barcharts
% Metrics: length density, branch density, fraction volume
regions = {'tiss','gyri','sulci','gm','wm','gm_sulci','wm_sulci',...
            'gm_gyri','wm_gyri'};
params = {'length_density','branch_density','fraction_volume','tortuosity'};
titles = {'Tissue Length Density','Tissue Branch Density',...
            'Tissue Volume Fraction'};
xlabels = {'Length Density (\mum^-^2)','Branch Density (\mum^-^3)',...
            'Volume Fraction (a.u.)'};
plot_names = {'length_density','branch_density','volume_fraction'};

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
        if j~=4
            plot_name = strcat(regions{ii},'_',plot_names{j});
            bar_chart(ad, cte, nc, subids, mpath, titles{j},'',...
                      xlabels{j}, plot_name);
            close;
        end
    end
end

%% Calculate ratio of sulci/gyri
% Cell array of parameters to compute ratio 
ratio_params = {'length_density','branch_density','fraction_volume'};
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
% Include the ratios in the regions
regions = {'tiss','gyri','sulci','gm','wm','gm_sulci','wm_sulci',...
            'gm_gyri','wm_gyri','sulci_gyri','wm_sulci_gyri','gm_sulci_gyri'};
% Calculate stats
pstats = metrics_stats(metrics, regions, params, groups, alpha, trend);

%% Generate table of p-values
% TODO: organize the p-values into a table. Call the function "make_ptable"

%{
% Pairwise comparison names
Pairs = {'AD vs CTE', 'AD vs. HC', 'CTE vs. HC'};
LengthDensity = ;
BranchDensity = ;
VolumeFraction = ;
Tortuosity = ;
%}

%% Generate histograms/violin plots of tortuosity
% TODO: generalize a function and call for all tortuosities

%% Box / Whisker Plots
%{
%%% total length (microns)
figure;
x1 = ad_total_len;
x2 = cte_total_len;
x3 = nc_total_len;
x = [x1; x2; x3];
g = [zeros(length(x1), 1); ones(length(x2), 1); ones(length(x3), 1).*2];
b = boxplot(x, g,'Labels',{'AD', 'CTE', 'NC'});
title('Total Length (\mum)')
set(gca, 'FontSize', 30)
set(b,'LineWidth',4)

%%% length density
figure;
x1 = ad_lenden;
x2 = cte_lenden;
x3 = nc_lenden;
x = [x1; x2; x3];
b = boxplot(x, g,'Labels',{'AD', 'CTE', 'NC'});
title({'Length Density','total vessel length / metric volume',''})
set(gca, 'FontSize', 30)
set(b,'LineWidth',4)

%%% length density (AD vs. NC)
figure;
x1 = ad_lenden;
x3 = nc_lenden;
x = [x1; x3];
gsub = [zeros(length(x1), 1); ones(length(x3), 1).*2];
b = boxplot(x, gsub,'Labels',{'AD', 'HC'});
title('Length Density');
ylabel('(10^-^7)(\mum^-^2)')
set(gca, 'FontSize', 80)
set(b,'LineWidth',4)

%%% total vessels
figure;
x1 = ad_nves;
x2 = cte_nves;
x3 = nc_nves;
x = [x1; x2; x3];
g = [zeros(length(x1), 1); ones(length(x2), 1); ones(length(x3), 1).*2];
b = boxplot(x, g,'Labels',{'AD', 'CTE', 'NC'});
title('Total Vessels')
set(gca, 'FontSize', 30)
set(b,'LineWidth',4)


%%% tortuosity (unitless)
figure;
x1 = ad_tort;
x2 = cte_tort;
x3 = nc_tort;
x = [x1; x2; x3];
g = [zeros(length(x1), 1); ones(length(x2), 1); ones(length(x3), 1).*2];
b = boxplot(x, g,'Labels',{'AD', 'CTE', 'NC'});
title('Tortuosity (\mum)')
set(gca, 'FontSize', 30)
set(b,'LineWidth',4)

%%% Branch Density (branches / mm^3)
figure;
x1 = ad_bden;
x2 = cte_bden;
x3 = nc_bden;
x = [x1; x2; x3];
g = [zeros(length(x1), 1); ones(length(x2), 1); ones(length(x3), 1).*2];
b = boxplot(x, g,'Labels',{'AD', 'CTE', 'NC'});
title('Branch Density (branches / mm^3)')
set(gca, 'FontSize', 30)
set(b,'LineWidth',4)

%%% Branch Density (CTE vs. NC)
figure;
x2 = cte_bden;
x3 = nc_bden;
x = [x2; x3];
gsub = [ones(length(x2), 1); ones(length(x3), 1).*2];
b = boxplot(x, gsub,'Labels',{'CTE', 'HC'});
title('Branch Density')
ylabel('(Branches / mm^3)')
set(gca, 'FontSize', 50)
set(b,'LineWidth',4)

%%% Fraction Volume (unitless)
figure;
x1 = ad_fvol;
x2 = cte_fvol;
x3 = nc_fvol;
x = [x1; x2; x3];
g = [zeros(length(x1), 1); ones(length(x2), 1); ones(length(x3), 1).*2];
b = boxplot(x, g,'Labels',{'AD', 'CTE', 'NC'});
title('Fraction Volume (Unitless)')
set(gca, 'FontSize', 30)
set(b,'LineWidth',4)

%%% Fraction Volume (AD vs. HC)
figure;
x1 = ad_fvol;
x3 = nc_fvol;
x = [x1; x3];
gsub = [zeros(length(x1), 1); ones(length(x3), 1).*2];
b = boxplot(x, gsub,'Labels',{'AD','HC'});
title('Fraction Volume')
ylabel('(a.u.)(10^-^3)')
yticklabels({'0.6', '1.2', '1.8'});
yticks([0.6, 1.2, 1.8] * 1e-3);
set(gca, 'FontSize', 80)
set(b,'LineWidth',4)

%}
