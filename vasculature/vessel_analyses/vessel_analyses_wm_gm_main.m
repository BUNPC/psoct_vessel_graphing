%% Vessel metrics - statistical analyses
%{
This script is for analyzing the length density in the dataset:
/projectnb/npbssmic/ns/Ann_Mckee_samples_10T

This script compares both the entire tissue volumes and the white vs gray
matter.

These steps must be performed before reanalyzing data:
- Recalculate metrics:
    - Diameter
    - Tortuosity
    - Length density (total length per volume)
    - Perivascular space (currently unable to differentiate)
    - Blood vessel branch density
        - # blood vessel branches / unit volume
        - # blood vessel branches / tree
    - Fraction Volume (volume of blood vessels / volume of image stack)
    - The sulcus to crest ratio for each of the following:
        - Gray matter branch density
        - Gray matter fraction volume
        - White matter branch density
        - White matter fraction volume
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

%% Data Directory
% Path to top-level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
% Volume and graph directories
voldir = '/dist_corrected/volume/ref/';
maskdir = '/dist_corrected/volume/ref/masks/';
graphdir = ['/dist_corrected/volume/combined_segs/' ...
            'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];
% Metrics output path
mpath = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
    'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];

%%% Subjects
ad_sub = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424'};
cte_sub = {'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_7126', 'CTE_8572'};
nc_sub = {'NC_6839','NC_6974','NC_8653','NC_21499','NC_301181'};
%%% Filenames
% Tissue mask filename (same for each subject)
t_mask = 'maskec.mat';
% White matter mask filename
wm_mask = 'wm_mask.mat';
% Gray matter mask filename
gm_mask = 'gm_mask.mat';
% graph matlab struct
graphname = 'seg_refined_masked_rmloop_graph_data.mat';
wm_graphname = 'seg_wm_refined_masked_rmloop_graph_data.mat';
gm_graphname = 'seg_gm_refined_masked_rmloop_graph_data.mat';

%%% Output filenames for metrics structs
% AD / CTE
ad_cte_fname = 'AD_CTE_metrics.mat';
ad_cte_fout = fullfile(mpath, ad_cte_fname);
ad_cte_wm_fout = fullfile(mpath, 'AD_CTE_wm_metrics.mat');
ad_cte_gm_fout = fullfile(mpath, 'AD_CTE_gm_metrics.mat');

% NC
nc_fname = 'NC_metrics.mat';
nc_fout = fullfile(mpath, nc_fname);
nc_wm_fout = fullfile(mpath, 'NC_wm_metrics.mat');
nc_gm_fout = fullfile(mpath, 'NC_gm_metrics.mat');

% ANOVA
anova_fname = 'ANOVA_ad_cte_nc.mat';
anova_fout = fullfile(mpath, anova_fname);

% Voxel dimensions (microns) and volume (cubic micron)
vox_dim = [12, 12, 15];
vox_vol = vox_dim(1) .* vox_dim(2) .* vox_dim(3);

%% Calculate metrics
% Metrics are: total length, length density, mean length, tortuosity

%%% Gray Matter
%{
% AD / CTE
subid = horzcat(ad_sub, cte_sub);
[met] = stats(dpath, subid, voldir, maskdir, graphdir, gm_graphname, gm_mask,...
                vox_vol, vox_dim);
save(ad_cte_gm_fout, 'met', '-v7.3');
% NC
[met] = stats(dpath, nc_sub, voldir, maskdir, graphdir, gm_graphname, gm_mask,...
                vox_vol, vox_dim);
save(nc_gm_fout, 'met', '-v7.3');
%}

%%% White Matter
%{
% AD / CTE
subid = horzcat(ad_sub, cte_sub);
[met] = stats(dpath, subid, voldir, maskdir, graphdir, wm_graphname, wm_mask,...
                vox_vol, vox_dim);
save(ad_cte_wm_fout, 'met', '-v7.3');
% NC
[met] = stats(dpath, nc_sub, voldir, maskdir, graphdir, wm_graphname, wm_mask,...
                vox_vol, vox_dim);
save(nc_wm_fout, 'met', '-v7.3');
%}

%%% Entire tissue volume
%{
% AD / CTE
subid = horzcat(ad_sub, cte_sub);
[met] = stats(dpath, subid, voldir, [], graphdir, graphname, t_mask,...
                vox_vol, vox_dim);
save(ad_cte_fout, 'met', '-v7.3');
% NC
[met] = stats(dpath, nc_sub, voldir, [], graphdir, graphname, t_mask,...
                vox_vol, vox_dim);
save(nc_fout, 'met', '-v7.3');
%}

%% Generate Barcharts of metrics for entire tissue volumes

%%% Filename prefix for WM and GM
% fprefix = 'volume';
% fprefix = 'wm';
fprefix = 'gm';

%%% Title prefix for WM and GM
% tprefix = 'Volume';
% tprefix = 'White Matter';
tprefix = 'Gray Matter';

%%% Indices in structure corresponding to AD and CTE
n_ad = length(ad_sub);
n_cte = length(cte_sub);
ad_idx = 1:n_ad;
cte_idx = (n_ad+1):(n_ad+n_cte);

%%% load metrics for NC (select either tissue, wm, or gm)
% nc = load(nc_fout);
% nc = load(nc_wm_fout);
nc = load(nc_gm_fout);
nc = nc.met;

%%% load metrics for AD & CTE (select either tissue, wm, or gm)
% met = load(ad_cte_fout);
% met = load(ad_cte_wm_fout);
met = load(ad_cte_gm_fout);
met = met.met;
% Separate CTE and AD
ad = met(ad_idx);
cte = met(cte_idx);

%%% total length (microns)
mout = fullfile(mpath, strcat(fprefix,'_total_length_bar'));
ad_total_length = vertcat(ad.total_length);
cte_total_length = vertcat(cte.total_length);
nc_total_length = vertcat(nc.total_length);
% Title
tstring = strcat(tprefix,' Total Vessel Length (\mum)');
% Plot barchart
barchart(ad_total_length, cte_total_length, nc_total_length,...
        ad_sub, cte_sub, nc_sub,...
        tstring,'Length (\mum)',mout);

%%% average length (microns)
mout = fullfile(mpath, strcat(fprefix,'_avg_length_bar'));
ad_avg_len = vertcat(ad.avg_length);
cte_avg_len = vertcat(cte.avg_length);
nc_avg_len = vertcat(nc.avg_length);
% Title
tstring = strcat(tprefix,' Average Vessel Length (\mum)');
% Plot barchart
barchart(ad_avg_len, cte_avg_len, nc_avg_len,...
        ad_sub, cte_sub, nc_sub,...
        tstring,'Length (\mum)',mout);

%%% length density
mout = fullfile(mpath, strcat(fprefix,'_length_density_bar'));
ad_lenden = vertcat(ad.length_density);
cte_lenden = vertcat(cte.length_density);
nc_lenden =  vertcat(nc.length_density);
% Title
tstring = strcat(tprefix,' Length Density');
% Plot barchart
barchart(ad_lenden, cte_lenden, nc_lenden,...
        ad_sub, cte_sub, nc_sub,...
        tstring,'Length/Volume (1/\mum^2)',mout);

%%% total vessels
mout = fullfile(mpath, strcat(fprefix,'_total_vessels_bar'));
ad_nves = vertcat(ad.total_vessels);
cte_nves = vertcat(cte.total_vessels);
nc_nves = vertcat(nc.total_vessels);
% Title
tstring = strcat(tprefix,' Total Vessels per Sample');
% Plot barchart
barchart(ad_nves, cte_nves, nc_nves,...
        ad_sub, cte_sub, nc_sub,...
        tstring,'Number of Vessels',mout);

%%% tortuosity (unitless)
mout = fullfile(mpath, strcat(fprefix,'_tortuosity_bar'));
ad_tort = vertcat(ad.tortuosity);
cte_tort = vertcat(cte.tortuosity);
nc_tort = vertcat(nc.tortuosity);
% Title
tstring = strcat(tprefix,' Tortuosity');
% Plot barchart
barchart(ad_tort, cte_tort, nc_tort,...
        ad_sub, cte_sub, nc_sub,...
        tstring,'a.u.',mout);

%%% Branch Density (branch / mm^3)
mout = fullfile(mpath, strcat(fprefix,'_branch_density_bar'));
% Create branch density arrays. Convert from branch/um^3 -> branch / mm^3
ad_bden = vertcat(ad.branchden) .* 1e9;
cte_bden = vertcat(cte.branchden) .* 1e9;
nc_bden = vertcat(nc.branchden) .* 1e9;
% Title
tstring = strcat(tprefix,' Branch Density (branches / mm^3)');
% Plot barchart
barchart(ad_bden, cte_bden, nc_bden,...
        ad_sub, cte_sub, nc_sub,...
        tstring,'Branches / mm^3',mout);

%%% Volume Fraction (unitless)
mout = fullfile(mpath, strcat(fprefix,'_volume_fraction_bar'));
% Create vessel fraction volume arrays.
ad_fvol = vertcat(ad.vol_frac);
cte_fvol = vertcat(cte.vol_frac);
nc_fvol = vertcat(nc.vol_frac);
% Title
tstring = strcat(tprefix,' Volume Fraction');
% Plot barchart
barchart(ad_fvol, cte_fvol, nc_fvol,...
        ad_sub, cte_sub, nc_sub,...
        tstring,'a.u.',mout);
%}

%% Box / Whisker Plots

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


%% Calculate ANOVA for: length density, total length, total # vessels

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











%% Function to calculate metrics
function [met] = stats(dpath, subid, voldir, maskdir, graphdir, graphname,...
                        t_mask, vox_vol, vox_dim)
% Initialize struct for storing metrics
met = struct();
for ii = 1:length(subid)
    %%% Load OCT tissue mask and calculate total metric volume
    % If the maskdir is empty, then load the entire tissue volume mask
    if isempty(maskdir)
        fullpath = fullfile(dpath, subid{ii}, voldir);
        filename = strcat(fullpath, t_mask);
        % Load the tissue mask
        tissue = load(filename);
        tissue = tissue.mask;
    else
        % Load the WM or GM mask
        fullpath = fullfile(dpath, subid{ii}, maskdir);
        filename = strcat(fullpath, t_mask);
        tissue = load(filename);
        try
            tissue = tissue.wm_mat;
        catch
            tissue = tissue.gm_mat;
        end
    end
    
    % Calculate total number of non-zero voxels in tissue sample
    tissue = logical(tissue);
    voxels = sum(tissue(:));
    % Convert voxels to metric volume (cubic microns)
    vol = voxels .* vox_vol;
    
    %%% Save subject ID to struct
    met(ii).subID = string(subid{ii});
    
    %%% Save total volume to struct
    met(ii).volume = vol;

    %%% Load Graph and Segmentation
    % Define entire filepath 
    fullpath = fullfile(dpath, subid{ii}, graphdir);
    filename = strcat(fullpath, graphname);
    %%% Load Data struct
    Data = load(filename, 'Data');
    % Load graph and seg
    graph = Data.Data.Graph;
    seg = Data.Data.angio;
    % Load end nodes
    endnodes = graph.segInfo.segEndNodes;

    %%% Calculate total length (microns) from graph
    len = graph.segInfo.segLen_um;
    len_tot = sum(len(:));
    met(ii).total_length = len_tot;
    
    %%% Calculate mean length (microns)
    len_avg = mean(len);
    met(ii).avg_length = len_avg;

    %%% Calculate length density
    len_density = len_tot ./ vol;
    met(ii).length_density = len_density;

    %%% Total number of vessels
    nves = length(len);
    met(ii).total_vessels = nves;

    %%% Segmentation volume (vessel volume)
    % Convert segmentation to logical
    seglog = logical(seg);
    % Take sum of logical segmentation (volume in voxels)
    segvol = sum(seglog(:));
    % Convert volume (voxels) to cubic microns
    segvol = segvol .* vox_vol;
    met(ii).segvol = segvol;

    %%% Volume Fraction of blood vessels
    met(ii).vol_frac = segvol ./ vol;

    %%% Branching density
    % Number of branches
    nb = Data.Data.Graph.nB;
    % Branch density (# branches / unit volume (cubic micron))
    branchden = sum(nb) / vol;
    met(ii).branchden = branchden;
    
    %%% tortuosity arc-chord ratio (curve length / euclidean)
    % Initalize matrix for storing tortuosity
    tort = zeros(nves, 1);
    for j=1:nves
        % convert segment end nodes to cartesian coordinate
        node1 = graph.nodes(endnodes(j,1), :);
        node2 = graph.nodes(endnodes(j,2), :);
        % Convert cartesian coordinate to offset in microns
        node1 = node1 .* vox_dim;
        node2 = node2 .* vox_dim;
        % Calcualte euclidean distance of segment
        euc = sqrt((node1(1) - node2(1)).^2 +...
                    (node1(2) - node2(2)).^2 +...
                    (node1(3) - node2(3)).^2);
        % Calculate tortuosity (single arc-chord ratio)
        tort(j) = len(j) ./ euc;
    end
    % Remove infinite tortuosity (loops)
    tort(tort==inf) = [];
    % save the entire tortuosity array
    met(ii).tort = tort;
    % Add to metrics structures
    met(ii).tortuosity = mean(tort);   
end
end

function barchart(ad_metric, cte_metric, nc_metric,...
                    ad_sub, cte_sub, nc_sub,...
                    ptitle, ytitle, figout)
% BARCHART create barchart of the metric
% INPUTS:
%   *_metric (double array): array of values for the metric
%   *_sub (cell array): cell of strings of the subject IDs
%   ptitle (string): title of figure
%   ytitle (string): y-axis label of figure
%   figout (string): directory/filename to save output

% Concatenate all groups together
metric = vertcat(ad_metric, cte_metric, nc_metric);

% Create the figure
figure('units','normalized','outerposition',[0 0 1 1])
x = categorical(horzcat(ad_sub, cte_sub, nc_sub));
b = bar(x, metric);
title(ptitle)
ylabel(ytitle)
xlabel('Subject ID')
set(gca, 'FontSize', 30)
% Set color of bars
b.FaceColor = 'flat';
b.CData(1:5,:) = [.5, 0, .5; .5, 0, .5; .5, 0, .5; .5, 0, .5; .5, 0, .5];
b.CData(end-4:end, :) =...
    [0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0];
% Save output
saveas(gca, figout, 'png')
close all;
end