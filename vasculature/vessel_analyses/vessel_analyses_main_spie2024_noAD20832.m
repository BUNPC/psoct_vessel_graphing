%% Vessel length density Analysis
%{
This script is for analyzing the length density in the dataset:
/projectnb/npbssmic/ns/Ann_Mckee_samples_10T

These steps must be performed before reanalyzing data:
- Select Frangi segmentation with best performance
    - Segment with voxel intensity threshold
    - Compare performance
- Apply loop removal to graph
- Segment the white and grey matter for each volume.
    - Apply threshold to the TIF file mus[#].tif located in:
        /projectnb/npbssmic/ns/Ann_Mckee_samples_55T/[subject_ID]/fitting/
    - This is the scattering map, where the number represents the
    scattering coefficient for a 150um slice. The same scattering
    coefficient is used for every axial voxel in the 150um slice.
    - Apply a threshold to this map to segment the white and grey matter.
    THen apply this mask to each layer in the 150um slice. Repeat this for
    every 150um in the entire tissue volume.
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
topdir = mydir(1:idcs(end));
addpath(genpath(topdir));

%% Data Directory on laptop

% Path to top-level directory
dpath = ['C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\' ...
    'test_data\Ann_Mckee_samples_10T\'];

% Metrics output path
mpath = ['C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\' ...
    'test_data\Ann_Mckee_samples_10T\metrics\spie_pho_west_2024\'];

% Volume directory + volume filename (same for each subject)
voldir = '\dist_corrected\volume\';
volname = 'ref_4ds_norm_inv.tif';

% Output filenames for metrics structs
ad_cte_fname = 'AD_CTE_metrics_noAD20832.mat';
ad_cte_fout = fullfile(mpath, ad_cte_fname);
nc_fname = 'NC_metrics.mat';
nc_fout = fullfile(mpath, nc_fname);
anova_fname = 'ANOVA_ad_cte_nc_noAD20832.mat';
anova_fout = fullfile(mpath, anova_fname);

%%% graph matlab struct
graphname = 'seg_maskedgraph_data.mat';
graphdir = '\combined_segs\gsigma_3-5-7_5-7-9_7-9-11\';

% Voxel dimensions (microns) and volume (cubic micron)
vox_dim = [12, 12, 15];
vox_vol = vox_dim(1) .* vox_dim(2) .* vox_dim(3);

%%% AD / CTE subject IDs and directories
% AD/CTE subject ID list for Ann_Mckee_samples_10T
% subid = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424',...
%          'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_7126'};
% subid = {'AD_10382', 'AD_20969', 'AD_21354', 'AD_21424',...
%          'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_7126'};

%%% NC subject IDs and directories
% Normal Control subject ID list for Ann_Mckee_samples_10T
subid = {'NC_301181','NC_21499','NC_6839','NC_6974','NC_7597','NC_8653'};


%% Calculate metrics (total length, length density, mean length, tortuosity)
%{
%%% Initialize struct for storing metrics
met = struct();

for ii = 1:length(subid)
    %%% Load psoct tissue volume and calculate total metric volume
    % Define entire filepath 
    fullpath = fullfile(dpath, subid{ii}, voldir);
    filename = strcat(fullpath, volname);
    % Convert .tif to .MAT
    tissue = TIFF2MAT(filename);
    % Invert so that non-tissue voxels are zeros
    tissue_inv = imcomplement(tissue);
    % Calculate total number of non-zero voxels in tissue sample
    tissue_logical = logical(tissue_inv);
    voxels = sum(tissue_logical(:));
    % Convert voxels to metric volume (cubic microns)
    vol = voxels .* vox_vol;
    
    %%% Save subject ID to struct
    met(ii).subID = string(subid{ii});
    
    %%% Save total volume to struct
    met(ii).volume = vol;

    %%% Load Graph and Segmentation
    % Define entire filepath 
    fullpath = fullfile(dpath, subid{ii}, voldir, graphdir);
    filename = strcat(fullpath, graphname);
    %%% Load Data struct
    Data = load(filename, 'Data');
    % Load graph and seg
    graph = Data.Data.Graph;
    seg = Data.Data.angio;
    % Load end nodes
    endnodes = graph.segInfo.segEndNodes;
    % Load position of all nodes
    nodepos = graph.segInfo.segPos;

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

% save(ad_cte_fout, 'met', '-v7.3');
save(nc_fout, 'met', '-v7.3');
%}

%% Generate Barcharts of metrics
%{
%%%% load metrics for AD & CTE
met = load(ad_cte_fout);
met = met.met;
% Indices in matrix corresponding to AD and CTE
ad_idx = 1:4;
cte_idx = 5:8;
% Separate CTE and AD
ad = met(ad_idx);
cte = met(cte_idx);

%%%% load metrics for NC
met = load(nc_fout);
nc = met.met;

%%% total length (microns)
ad_total_len = vertcat(ad.total_length);
cte_total_len = vertcat(cte.total_length);
nc_total_len = vertcat(nc.total_length);
total_len = vertcat(ad_total_len, cte_total_len, nc_total_len);
figure('units','normalized','outerposition',[0 0 1 1])
x = categorical(...
    {'AD 10382', 'AD 20969', 'AD 21354','AD 21424',...
    'CTE 6489', 'CTE 6912', 'CTE 7019', 'CTE 7126',...
    'NC 301181','NC 21499','NC 6839','NC 6974','NC 7597','NC 8653'});
b = bar(x, total_len);
title('Total Vessel Length (\mum)')
ylabel('Length (\mum)')
xlabel('Subject ID')
set(gca, 'FontSize', 30)
% Set color of bars
b.FaceColor = 'flat';
b.CData(1:4,:) = [.5, 0, .5; .5, 0, .5; .5, 0, .5; .5, 0, .5];
b.CData(end-5:end, :) =...
    [0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0;];
% Save output
mout = fullfile(mpath, 'total_length_bar_noAD20832');
saveas(gca, mout, 'png')


%%% average length (microns)
ad_avg_len = vertcat(ad.avg_length);
cte_avg_len = vertcat(cte.avg_length);
nc_avg_len = vertcat(nc.avg_length);
avg_len = vertcat(ad_avg_len, cte_avg_len, nc_avg_len);
figure('units','normalized','outerposition',[0 0 1 1])
x = categorical(...
    {'AD 10382', 'AD 20969', 'AD 21354','AD 21424',...
    'CTE 6489', 'CTE 6912', 'CTE 7019', 'CTE 7126',...
    'NC 301181','NC 21499','NC 6839','NC 6974','NC 7597','NC 8653'});
b = bar(x, avg_len);
title('Average Vessel Length (\mum)')
ylabel('Length (\mum)')
xlabel('Subject ID')
set(gca, 'FontSize', 30)
% Set color of bars
b.FaceColor = 'flat';
b.CData(1:4,:) = [.5, 0, .5; .5, 0, .5; .5, 0, .5; .5, 0, .5];
b.CData(end-5:end, :) =...
    [0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0;];
% Save output
mout = fullfile(mpath, 'avg_length_bar_noAD20832');
saveas(gca, mout, 'png')


%%% length density
ad_lenden = vertcat(ad.length_density);
cte_lenden = vertcat(cte.length_density);
nc_lenden =  vertcat(nc.length_density);
lenden = vertcat(ad_lenden, cte_lenden, nc_lenden);
figure('units','normalized','outerposition',[0 0 1 1])
x = categorical(...
    {'AD 10382', 'AD 20969', 'AD 21354','AD 21424',...
    'CTE 6489', 'CTE 6912', 'CTE 7019', 'CTE 7126',...
    'NC 301181','NC 21499','NC 6839','NC 6974','NC 7597','NC 8653'});
b = bar(x, lenden);
title('Vessel Length Density (\mum / \mu^3)')
ylabel('Length Density (\mum^2)')
xlabel('Subject ID')
set(gca, 'FontSize', 30)
% Set color of bars
b.FaceColor = 'flat';
b.CData(1:4,:) = [.5, 0, .5; .5, 0, .5; .5, 0, .5; .5, 0, .5];
b.CData(end-5:end, :) =...
    [0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0;];
% Save output
mout = fullfile(mpath, 'length_density_bar_noAD20832');
saveas(gca, mout, 'png')

%%% total vessels
ad_nves = vertcat(ad.total_vessels);
cte_nves = vertcat(cte.total_vessels);
nc_nves = vertcat(nc.total_vessels);
nves = vertcat(ad_nves, cte_nves, nc_nves);
figure('units','normalized','outerposition',[0 0 1 1])
x = categorical(...
    {'AD 10382', 'AD 20969', 'AD 21354','AD 21424',...
    'CTE 6489', 'CTE 6912', 'CTE 7019', 'CTE 7126',...
    'NC 301181','NC 21499','NC 6839','NC 6974','NC 7597','NC 8653'});
b = bar(x, nves);
title('Total Vessels per Sample')
ylabel('Number of Vessels')
xlabel('Subject ID')
set(gca, 'FontSize', 30)
% Set color of bars
b.FaceColor = 'flat';
b.CData(1:4,:) = [.5, 0, .5; .5, 0, .5; .5, 0, .5; .5, 0, .5];
b.CData(end-5:end, :) =...
    [0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0;];
% Save output
mout = fullfile(mpath, 'total_vessels_bar_noAD20832');
saveas(gca, mout, 'png')

%%% tortuosity (unitless)
%{
ad_tort = vertcat(ad.tortuosity);
cte_tort = vertcat(cte.tortuosity);
nc_tort = vertcat(nc.tortuosity);
tort = vertcat(ad_tort, cte_tort, nc_tort );
figure('units','normalized','outerposition',[0 0 1 1])
x = categorical(...
    {'AD 10382', 'AD 20969', 'AD 21354','AD 21424',...
    'CTE 6489', 'CTE 6912', 'CTE 7019', 'CTE 7126',...
    'NC 301181','NC 21499','NC 6839','NC 6974','NC 7597','NC 8653'});
b = bar(x, tort);
title('Tortuosity (curve length / euclidean distance)')
ylabel('Tortuosity (arc-chord ratio)')
xlabel('Subject ID')
set(gca, 'FontSize', 30)
% Set color of bars
b.FaceColor = 'flat';
b.CData(1:4,:) = [.5, 0, .5; .5, 0, .5; .5, 0, .5; .5, 0, .5];
b.CData(end-5:end, :) =...
    [0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0;];
% Save output
mout = fullfile(mpath, 'tortuosity_bar_noAD20832');
saveas(gca, mout, 'png')
%}

%%% Branch Density (branch / mm^3)
% Create branch density arrays. Convert from branch/um^3 -> branch / mm^3
ad_bden = vertcat(ad.branchden) .* 1e9;
cte_bden = vertcat(cte.branchden) .* 1e9;
nc_bden = vertcat(nc.branchden) .* 1e9;
% Create bar chart
bden = vertcat(ad_bden, cte_bden, nc_bden);
figure('units','normalized','outerposition',[0 0 1 1])
x = categorical(...
    {'AD 10382', 'AD 20969', 'AD 21354','AD 21424',...
    'CTE 6489', 'CTE 6912', 'CTE 7019', 'CTE 7126',...
    'NC 301181','NC 21499','NC 6839','NC 6974','NC 7597','NC 8653'});
b = bar(x, bden);
title('Branch Density (branches / mm^3)')
ylabel('Branches / mm^3')
xlabel('Subject ID')
set(gca, 'FontSize', 30)
% Set color of bars
b.FaceColor = 'flat';
b.CData(1:4,:) = [.5, 0, .5; .5, 0, .5; .5, 0, .5; .5, 0, .5];
b.CData(end-5:end, :) =...
    [0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0;];
% Save output
mout = fullfile(mpath, 'branch_density_bar_noAD20832');
saveas(gca, mout, 'png')

%%% Fraction Volume (unitless)
% Create vessel fraction volume arrays.
ad_fvol = vertcat(ad.segvol) ./ vertcat(ad.volume);
cte_fvol = vertcat(cte.segvol) ./ vertcat(cte.volume);
nc_fvol = vertcat(nc.segvol) ./ vertcat(nc.volume);
% Create bar chart
fvol = vertcat(ad_fvol, cte_fvol, nc_fvol);
figure('units','normalized','outerposition',[0 0 1 1])
x = categorical(...
    {'AD 10382', 'AD 20969', 'AD 21354','AD 21424',...
    'CTE 6489', 'CTE 6912', 'CTE 7019', 'CTE 7126',...
    'NC 301181','NC 21499','NC 6839','NC 6974','NC 7597','NC 8653'});
b = bar(x, fvol);
title('Vessel Fraction Volume')
ylabel('um^3 / um^3')
xlabel('Subject ID')
set(gca, 'FontSize', 30)
% Set color of bars
b.FaceColor = 'flat';
b.CData(1:4,:) = [.5, 0, .5; .5, 0, .5; .5, 0, .5; .5, 0, .5];
b.CData(end-5:end, :) =...
    [0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0;];
% Save output
mout = fullfile(mpath, 'vessel_fraction_volume_bar_noAD20832');
saveas(gca, mout, 'png')
%}

%% Calcuate average + std. dev of metrics
%{
%%%% load metrics for AD & CTE
met = load(ad_cte_fout);
met = met.met;
% Indices in matrix corresponding to AD and CTE
ad_idx = 1:4;
cte_idx = 5:8;
% Separate CTE and AD
ad = met(ad_idx);
cte = met(cte_idx);

%%%% load metrics for NC
met = load(nc_fout);
nc = met.met;

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

%}


%% Calculate ANOVA for: length density, total length, total # vessels

%%% Load the AD,CTE struct and the NC struct
% load metrics for AD & CTE
met = load(ad_cte_fout);
met = met.met;
% Indices in matrix corresponding to AD and CTE
ad_idx = 1:4;
cte_idx = 5:8;
% Separate CTE and AD
ad = met(ad_idx);
cte = met(cte_idx);
% load metrics for NC ('NC_6839','NC_8095', 'NC_8653', 'NC_21499')
met = load(nc_fout);
nc = met.met;

%%% Define sample size for each group
n_ad = 4;
n_cte = 4;
n_nc = 6;

%%% Define arrays for unbalanced ANOVA
% Factor name arrays
g_ad_nc = [repmat("AD",n_ad,1); repmat("NC",n_nc,1)];
g_cte_nc = [repmat("CTE",n_cte,1); repmat("NC",n_nc,1)];
g_ad_cte = [repmat("AD",n_ad,1); repmat("CTE",n_cte,1)];

% Tortuosity
ad_tort = vertcat(ad.tort);
nc_tort = vertcat(nc.tort);
cte_tort = vertcat(cte.tort);
ad_nc_tort = [ad_tort, nc_tort];
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

% Total length
aov.ad_nc_lentot = anova1(ad_nc_len, g_ad_nc);
title('AD vs NC Length Total')
aov.cte_nc_lentot = anova1(cte_nc_len, g_cte_nc);
title('CTE vs NC Length Total')
aov.ad_cte_lentot = anova1(ad_cte_len, g_ad_cte);
title('AD vs CTE Length Total')

% Total number vessels
aov.ad_nc_nves = anova1(ad_nc_nves, g_ad_nc);
title('AD vs NC # Vessels')
aov.cte_nc_nves = anova1(cte_nc_nves, g_cte_nc);
title('CTE vs NC # Vessels')
aov.ad_cte_nves = anova1(ad_cte_nves, g_ad_cte);
title('AD vs CTE # Vessels')

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
age_nc = [59;88;71;73;69;80];
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












