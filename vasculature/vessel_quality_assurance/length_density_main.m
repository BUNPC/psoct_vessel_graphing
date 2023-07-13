%% Vessel length density Analysis
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
topdir = mydir(1:idcs(end));
addpath(genpath(topdir));

%% Data Directory on laptop

% Path to top-level directory
dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\';

% Metrics output path
mpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\metrics\';

% Volume directory + volume filename (same for each subject)
voldir = '\dist_corrected\volume\';
volname = 'ref_4ds_norm_inv.tif';

% Output filenames for metrics structs
ad_cte_fname = 'AD_CTE_metrics.mat';
ad_cte_fout = fullfile(dpath, ad_cte_fname);
nc_fname = 'NC_metrics.mat';
nc_fout = fullfile(dpath, nc_fname);

% Name of graph matlab struct
graphname = 'ref_4ds_norm_inv_segment_pmin_0.23_mask_40_graph_data.mat';

% Voxel dimensions (microns) and volume (cubic micron)
vox_dim = [12, 12, 15];
vox_vol = vox_dim(1) .* vox_dim(2) .* vox_dim(3);

%%% AD / CTE subject IDs and directories
% AD/CTE subject ID list for Ann_Mckee_samples_10T
subid = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_7126'};
% AD/CTE segment directory + filename
segdir = {'gsigma_1-2-3-4-5_gsize_5--9-13-17-21\',...
    'gsigma_1-2-3-4-5_gsize_5--9-13-17-21\',...
    'gsigma_1-2-3-4-5_gsize_5--9-13-17-21\',...
    'gsigma_1-2-3-4-5_gsize_5--9-13-17-21\',...
    'gsigma_1-2-3-4-5_gsize_5--9-13-17-21\',...
    'gsigma_1-2-3-4-5_gsize_5--9-13-17-21\',...
    'gsigma_1-2-3-4-5_gsize_5--9-13-17-21\',...
    'gsigma_1-2-3-4-5_gsize_5--9-13-17-21\',...
    'gsigma_3-5-7-9_gsize_13-21-29-37\'};

%%% NC subject IDs and directories
%{
% Normal Control subject ID list for Ann_Mckee_samples_10T
subid = {'NC_6047', 'NC_6839','NC_8095', 'NC_8653', 'NC_21499'};
% AD/CTE segment directory + filename
segdir = {
    'gsigma_1-2-3-4-5_gsize_5--9-13-17-21\',...
    'gsigma_1-2-3-4-5_gsize_5--9-13-17-21\',...
    'gsigma_1-3-5_gsize_5-13-21\',...
    'gsigma_1-3-5_gsize_5-13-21\',...
    'gsigma_1-2-3-4-5_gsize_5--9-13-17-21\'};
%}


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
    % Save total volume to struct
    met(ii).volume = vol;

    %%% Load segmentation
    % Define entire filepath 
    fullpath = fullfile(dpath, subid{ii}, voldir, segdir{ii});
    filename = strcat(fullpath, graphname);
    %%% Load Data struct
    Data = load(filename, 'Data');
    % Load grapha and seg
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
    
    %%% tortuosity arc-chord ratio (curve length / euclidean)
    % Initalize matrix for storing tortuosity
    tort = zeros(nves, 1);
    for j=1:nves
        % convert segment end nodes to cartesian coordinate
        node1 = graph.nodes(endnodes(j,1), :);
        node2 = graph.nodes(endnodes(j,2), :);
        % Calcualte euclidean distance of segment
        euc = sqrt((node1(1) - node2(1)).^2 +...
                    (node1(2) - node2(2)).^2 +...
                    (node1(3) - node2(3)).^2);
        % Calculate tortuosity (single arc-chord ratio)
        tort(j) = len(j) ./ euc;
    end
    % Add to metrics structures
    met(ii).tortuosity = mean(tort);   
end

save(ad_cte_fout, 'met', '-v7.3');
% save(nc_fout, 'met', '-v7.3');
%}

%% Generate Barcharts of metrics

%%%% load metrics for AD & CTE
met = load(ad_cte_fout);
met = met.met;
% Indices in matrix corresponding to AD and CTE
ad_idx = 1:5;
cte_idx = 6:9;
% Separate CTE and AD
ad = met(ad_idx);
cte = met(cte_idx);

%%%% load metrics for NC
met = load(nc_fout);
met = met.met;
nc = met(2:end);


%%% total length (microns)
ad_total_len = vertcat(ad.total_length);
cte_total_len = vertcat(cte.total_length);
nc_total_len = vertcat(nc.total_length);
total_len = vertcat(ad_total_len, cte_total_len, nc_total_len);
figure('units','normalized','outerposition',[0 0 1 1])
x = categorical({'AD 10382', 'AD 20832', 'AD 20969', 'AD 21354',...
                'AD 21424','CTE 6489', 'CTE 6912', 'CTE 7019', 'CTE 7126',...
                'NC 6839','NC 8095', 'NC 8653', 'NC 21499'});
b = bar(x, total_len);
title('Total Vasculature Length (\mum)')
ylabel('Vasculature Total Length (\mum)')
xlabel('Subject ID')
set(gca, 'FontSize', 30)
% Set color of bars
b.FaceColor = 'flat';
b.CData(1:4,:) = [.5, 0, .5; .5, 0, .5; .5, 0, .5; .5, 0, .5];
b.CData(end-3:end, :) = [0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0];
% Save output
mout = fullfile(mpath, 'total_length_bar');
saveas(gca, mout, 'png')


%%% average length (microns)
ad_avg_len = vertcat(ad.avg_length);
cte_avg_len = vertcat(cte.avg_length);
nc_avg_len = vertcat(nc.avg_length);
avg_len = vertcat(ad_avg_len, cte_avg_len, nc_avg_len);
figure('units','normalized','outerposition',[0 0 1 1])
x = categorical({'AD 10382', 'AD 20832', 'AD 20969', 'AD 21354',...
                'AD 21424','CTE 6489', 'CTE 6912', 'CTE 7019', 'CTE 7126',...
                'NC 6839','NC 8095', 'NC 8653', 'NC 21499'});
b = bar(x, avg_len);
title('Average Vessel Length (\mum)')
ylabel('Vasculature Average Length (\mum)')
xlabel('Subject ID')
set(gca, 'FontSize', 30)
% Set color of bars
b.FaceColor = 'flat';
b.CData(1:4,:) = [.5, 0, .5; .5, 0, .5; .5, 0, .5; .5, 0, .5];
b.CData(end-3:end, :) = [0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0];
% Save output
mout = fullfile(mpath, 'avg_length_bar');
saveas(gca, mout, 'png')


%%% length density
ad_lenden = vertcat(ad.length_density);
cte_lenden = vertcat(cte.length_density);
nc_lenden =  vertcat(nc.length_density);
lenden = vertcat(ad_lenden, cte_lenden, nc_lenden);
figure('units','normalized','outerposition',[0 0 1 1])
x = categorical({'AD 10382', 'AD 20832', 'AD 20969', 'AD 21354',...
                'AD 21424','CTE 6489', 'CTE 6912', 'CTE 7019', 'CTE 7126',...
                'NC 6839','NC 8095', 'NC 8653', 'NC 21499'});
b = bar(x, lenden);
title('Length Density (\mum / \mu^3)')
ylabel('Vasculature Length Density (\mum^2)')
xlabel('Subject ID')
set(gca, 'FontSize', 30)
% Set color of bars
b.FaceColor = 'flat';
b.CData(1:4,:) = [.5, 0, .5; .5, 0, .5; .5, 0, .5; .5, 0, .5];
b.CData(end-3:end, :) = [0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0];
% Save output
mout = fullfile(mpath, 'length_density_bar');
saveas(gca, mout, 'png')

%%% total vessels
ad_nves = vertcat(ad.total_vessels);
cte_nves = vertcat(cte.total_vessels);
nc_nves = vertcat(nc.total_vessels);
nves = vertcat(ad_nves, cte_nves, nc_nves);
figure('units','normalized','outerposition',[0 0 1 1])
x = categorical({'AD 10382', 'AD 20832', 'AD 20969', 'AD 21354',...
                'AD 21424','CTE 6489', 'CTE 6912', 'CTE 7019', 'CTE 7126',...
                'NC 6839','NC 8095', 'NC 8653', 'NC 21499'});
b = bar(x, nves);
title('Total # Vessels per Sample')
ylabel('Number of Vessels')
xlabel('Subject ID')
set(gca, 'FontSize', 30)
% Set color of bars
b.FaceColor = 'flat';
b.CData(1:4,:) = [.5, 0, .5; .5, 0, .5; .5, 0, .5; .5, 0, .5];
b.CData(end-3:end, :) = [0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0];
% Save output
mout = fullfile(mpath, 'total_vessels_bar');
saveas(gca, mout, 'png')
%}

%%% tortuosity (unitless)
ad_tort = vertcat(ad.tortuosity);
cte_tort = vertcat(cte.tortuosity);
nc_tort = vertcat(nc.tortuosity);
tort = vertcat(ad_tort, cte_tort, nc_tort );
figure('units','normalized','outerposition',[0 0 1 1])
x = categorical({'AD 10382', 'AD 20832', 'AD 20969', 'AD 21354',...
                'AD 21424','CTE 6489', 'CTE 6912', 'CTE 7019', 'CTE 7126',...
                'NC 6839','NC 8095', 'NC 8653', 'NC 21499'});
b = bar(x, tort);
title('Tortuosity (curve length / euclidean distance)')
ylabel('Tortuosity (arc-chord ratio)')
xlabel('Subject ID')
set(gca, 'FontSize', 30)
% Set color of bars
b.FaceColor = 'flat';
b.CData(1:4,:) = [.5, 0, .5; .5, 0, .5; .5, 0, .5; .5, 0, .5];
b.CData(end-3:end, :) = [0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0; 0, 0.5, 0];
% Save output
mout = fullfile(mpath, 'tortuosity_bar');
saveas(gca, mout, 'png')


%% Calcuate average + std. dev of metrics

%%%% load metrics for AD & CTE
met = load(ad_cte_fout);
met = met.met;
% Indices in matrix corresponding to AD and CTE
ad_idx = 1:5;
cte_idx = 6:9;
% Separate CTE and AD
ad = met(ad_idx);
cte = met(cte_idx);

%%%% load metrics for NC
met = load(nc_fout);
met = met.met;
nc = met(2:end);

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












