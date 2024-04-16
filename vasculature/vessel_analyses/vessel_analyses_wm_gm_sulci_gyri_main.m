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

%% Barcharts
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

%% Statistical Hypothesis Testing
% Trend threshold
th_trend = 0.10;
% Significant Difference threshold
th_sigdif = 0.05;
% Calculate stats
pstats = metrics_stats(metrics, regions, params, th_trend, th_sigdif);

%% Calculate ANOVA for: length density, branch density, fraction volume

%%% Load the AD,CTE struct and the NC struct
% load metrics for AD & CTE
met = load(ad_cte_fout);
met = met.met;
% Separate CTE and AD
ad = met(ad_idx);
cte = met(cte_idx);
% load metrics for NC ('NC_6839','NC_8095', 'NC_8653', 'NC_21499')
met = load(nc_fout);
nc = met.met;

%%% Define sample size for each group
n_ad = 5;
n_cte = 4;
n_nc = 4;

%%% Define arrays for unbalanced ANOVA
% Factor name arrays
g_ad_nc = [repmat("AD",n_ad,1); repmat("NC",n_nc,1)];
g_cte_nc = [repmat("CTE",n_cte,1); repmat("NC",n_nc,1)];
g_ad_cte = [repmat("AD",n_ad,1); repmat("CTE",n_cte,1)];

% Tortuosity
ad_tort = vertcat(ad.tort);
nc_tort = vertcat(nc.tort);
cte_tort = vertcat(cte.tort);
ad_nc_tort = [ad_tort', nc_tort'];
cte_nc_tort = [cte_tort', nc_tort'];
ad_cte_tort = [ad_tort',cte_tort'];
% Group labels for unbalanced
g_ad_nc_tort = [repmat("AD",length(ad_tort),1); repmat("NC",length(nc_tort),1)];
g_cte_nc_tort = [repmat("CTE",length(cte_tort),1); repmat("NC",length(nc_tort),1)];
g_ad_cte_tort = [repmat("AD",length(ad_tort),1); repmat("CTE",length(cte_tort),1)];

% Length density arrays
ad_nc_den = [ad_lenden', nc_lenden'];
cte_nc_den = [cte_lenden', nc_lenden'];
ad_cte_den = [ad_lenden', cte_lenden'];

% Total length arrays
ad_nc_len = [ad_total_len', nc_total_len'];
cte_nc_len = [cte_total_len', nc_total_len'];
ad_cte_len = [ad_total_len', cte_total_len'];

% Total number vessels
ad_nc_nves = [ad_nves', nc_nves'];
cte_nc_nves = [cte_nves', nc_nves'];
ad_cte_nves = [ad_nves',cte_nves'];

% Branch Density
ad_nc_bden = [ad_bden', nc_bden'];
cte_nc_bden = [cte_bden', nc_bden'];
ad_cte_bden = [ad_bden',cte_bden'];

% Fraction Volume
ad_nc_fvol = [ad_fvol', nc_fvol'];
cte_nc_fvol = [cte_fvol', nc_fvol'];
ad_cte_fvol = [ad_fvol',cte_fvol'];

%%% Perform one-way unbalanced ANOVA (AD vs. NC, CTE vs. NC)
% Tortuosity
aov.ad_nc_tort = anova1(ad_nc_tort, g_ad_nc_tort);
title('AD vs NC Tortuosity')
aov.cte_nc_tort = anova1(cte_nc_tort, g_cte_nc_tort);
title('CTE vs NC Tortuosity')
aov.ad_cte_tort = anova1(ad_cte_tort, g_ad_cte_tort);
title('AD vs CTE Tortuosity')

% Length density
aov.ad_nc_lenden = anova1(ad_nc_den, g_ad_nc);
title('AD vs NC Length Density')
aov.cte_nc_lenden = anova1(cte_nc_den, g_cte_nc);
title('CTE vs NC Length Density')
aov.ad_cte_lenden = anova1(ad_cte_den, g_ad_cte);
title('AD vs CTE Length Density')

% Branch Density
aov.ad_nc_bden = anova1(ad_nc_bden, g_ad_nc);
title('AD vs NC Branch Density')
aov.cte_nc_bden = anova1(cte_nc_bden, g_cte_nc);
title('CTE vs NC Branch Density')
aov.ad_cte_bden = anova1(ad_cte_bden, g_ad_cte);
title('AD vs CTE Branch Density')

% Fraction Volume
aov.ad_nc_fvol = anova1(ad_nc_fvol, g_ad_nc);
title('AD vs NC Fraction Volume')
aov.cte_nc_fvol = anova1(cte_nc_fvol, g_cte_nc);
title('CTE vs NC Fraction Volume')
aov.ad_cte_fvol = anova1(ad_cte_fvol, g_ad_cte);
title('AD vs CTE Fraction Volume')

%%% Save the ANOVA
save(anova_fout, 'aov', '-v7.3');


%% Calculate average + std of age

%{
AD subjects = AD_10382, AD_20969, AD_21354, AD_21424
CTE subjects = CTE_6489, CTE_6912, CTE_7019, CTE_7126
NC subjects = NC_301181, NC_21499, NC_6839, NC_6974, NC_7597, NC_8653
%}

age_ad = [84;76;86;83];
age_cte = [81; 78; 75; 86];
age_nc = [59;88;71;69;80];
age_all_groups = [age_ad; age_cte; age_nc];

% AD
age.ad_mean = mean(age_ad);
age.ad_std = std(age_ad);

% CTE
age.cte_mean = mean(age_cte);
age.cte_std = std(age_cte);

% NC
age.nc_mean = mean(age_nc);
age.nc_std = std(age_nc);

% Overall age stats
age.mean = mean(age_all_groups);
age.std = std(age_all_groups);
age.min = min(age_all_groups);
age.max = max(age_all_groups);
%}