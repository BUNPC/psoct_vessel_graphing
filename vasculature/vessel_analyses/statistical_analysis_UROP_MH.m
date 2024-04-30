%% Main script for kmeans, PCA, and graphical data analyses on CAA cases
clear;
clc;
close all;

%% Managing environment and paths
%{
% Add path to wokring directory and .mat data struct
addpath('/Users/jamieedelists/Desktop/spring_UROP_analyses');
addpath('/Users/jamieedelists/Desktop/spring_UROP_analyses/violin/');

% Load the data struct
struct_subs = load('/Users/jamieedelists/Desktop/spring_UROP_analyses/caa_loops_metrics.mat');
%}
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

% Directory to save outputs
output_dir = '/projectnb/npbssmic/ns/CAA/metrics/crop300_rmloop';

% Load the data struct
struct_subs = load('/projectnb/npbssmic/ns/CAA/metrics/caa_crop300_rmloops_metrics.mat');

%}
%% Sorting each metrics into it's own array

% Names of the fields corresponding to subjects in the struct for easier indexing
subjects = fieldnames(struct_subs.metrics);

% Initialize/preallocate arrays to store metrics for each subject
len_den = zeros(1, numel(subjects));
branch_den = zeros(1, numel(subjects));
frac_vol = zeros(1, numel(subjects));
tort_mean = zeros(1, numel(subjects));
tort_med = zeros(1, numel(subjects));
tort_mode = zeros(1, numel(subjects));

% Placing metrics into seperate arrays for concatenation (except tort_array)
% Extract metrics for each subject
for i = 1:numel(subjects)
    len_den(i) = struct_subs.metrics.(subjects{i}).length_density;
    branch_den(i) = struct_subs.metrics.(subjects{i}).branch_density;
    frac_vol(i) = struct_subs.metrics.(subjects{i}).fraction_volume;
    tort_mean(i) = struct_subs.metrics.(subjects{i}).tortuosity_mean;
    tort_med(i) = struct_subs.metrics.(subjects{i}).tortuosity_med;
    tort_mode(i) = struct_subs.metrics.(subjects{i}).tortuosity_mode;
end
%% Calling kMeans on all metrics (no tort array)

k = 3;  % CHANGE as needed (number of clusters)

% kMeans matrix variables and transpose for colunms to represent variables and rows for different cases
kMeans_X = [len_den; branch_den; frac_vol; tort_mean; tort_med; tort_mode]';

% rows should be 8 (number of subjects) and cols should be 3 + 3(from tort)
[idx, C] = kmeans(kMeans_X, k);

%% Create box/whisker plots kmeans1 (without PCA)
% Combine metrics into a matrix. Rows = metrics. Columns = samples.
metrics_mat = vertcat(len_den, branch_den, frac_vol, tort_mean,...
                    tort_med, tort_mode);
% Create box/whisker plot
fstring = 'kmeans1_';
group_box_whisker(metrics_mat, idx, output_dir, fstring)


%% Extracting principle components (running PCA)

% Perform PCA on the data
[coeff, score, latent, tsquared, explained] = pca(kMeans_X);

% Choose the number of components to keep based on explained variance
% Find components that explain at least 85% of the variance
explained_sum = cumsum(explained);
% Find components that explain at least (threshold)% of the variance
threshold = 85;
numComponents = find(explained_sum >= threshold, 1, 'first');

% Reduce data dimensions based on the number of principal componets (first
% row or first 2 rows)
reducedData = score(:, 1:numComponents);

% Perform k-means clustering on the reduced data
k = 3; % Number of clusters
[idx2, C2] = kmeans(reducedData, k);

%% Create box/whisker plot with the second k-Means clustering
fstring = 'kmeans2_';
group_box_whisker(metrics_mat, idx2, output_dir,fstring)

%% Graphing the kMeans clusters of principal components (NEEDS FIXING WITH GRAPH)
% Display the results
disp('Cluster Indices:');
disp(idx2);

% Check the number of principal components and visualize accordingly
if numComponents == 1
    figure;
    scatter(reducedData(:,1), ones(size(reducedData, 1), 1), 10, idx, 'filled');
    xlabel('Principal Component 1');
    title('PCA followed by k-Means Clustering (1D)');
elseif numComponents == 2
    figure;
    gscatter(reducedData(:,1), reducedData(:,2), idx);
    xlabel('Principal Component 1');
    ylabel('Principal Component 2');
    title('PCA followed by k-Means Clustering (2D)');
elseif numComponents == 3
    figure;
    scatter3(reducedData(:,1), reducedData(:,2), reducedData(:,3), 10, idx, 'filled');
    xlabel('Principal Component 1');
    ylabel('Principal Component 2');
    zlabel('Principal Component 3');
    view(3);
    title('PCA followed by k-Means Clustering (3D)');
end

%% Box and Whiskers plot for length_den, branch_den, frac_vol

%Ranges between the different metrics are too large, not good for
%visualizing together
%box_whiskers_data_1 = [(len_den'), (branch_den')];

figure;

% Subplot for length_density
subplot(1,3,1); 
boxplot(len_den', 'Labels', {'Length Density'});
title('Length Density');
ylabel('Length Density');

% Subplot for branch_density
subplot(1,3,2); 
boxplot(branch_den', 'Labels', {'Branch Density'});
title('Branch Density');
ylabel('Branch Density');

% Subplot for fraction_volume
subplot(1,3,3); 
boxplot(frac_vol', 'Labels', {'Fraction Volume'});
title('Fraction Volume');
ylabel('Fraction Volume');

% Main title for entire figure
%sgtitle('Box and Whiskers Plots for Different Metrics'); 

%% Violin plots for tortuosity arrays

% Find the length of longest tortuosity array
tort_lengths = zeros(numel(subjects), 1);

for i = 1:numel(subjects)
    tort_lengths(i) = numel(struct_subs.metrics.(subjects{i}).tortuosity_array);
end

max_L = max(tort_lengths);

% Preallocate Y matrix for violin plots of columns (tortuosity arrays)
Y = NaN(max_L, numel(subjects));

% Populate tortuosity matrix: each column is tort array for a given subject
for i = 1:numel(subjects)
    currentLength = numel(struct_subs.metrics.(subjects{i}).tortuosity_array);
    Y(1:currentLength, i) = struct_subs.metrics.(subjects{i}).tortuosity_array;
end

figure; 

%Violin plot with distribution outlined:
[h, L, MX, MED, bw] = violin(Y, 'facecolor', [0.5, 0.5, 0.5], 'edgecolor', 'k', 'mc', 'k', 'medc', 'r'); 

%Violin plot with distribution not outlined:
%[h, L, MX, MED, bw] = violin(Y, 'facecolor', [0.5, 0.5, 0.5], 'edgecolor', 'none', 'mc', 'k', 'medc', 'r');

% Customizing windowing
ylim([0.8, 3.5]);  % REPLACE lower and upper bounds as needed

% Plotting violins
xlabel('Subjects');
ylabel('Tortuosity');
title('Violin Plots of Tortuosity For Each Subject');
% xlabels (x-ticks) set to the subject names
set(gca, 'XTick', 1:numel(subjects), 'XTickLabel', subjects);

%% Function to create a box/whisker plot

function group_box_whisker(metrics_mat, idx, output_dir, fstring)
% Generate a box/whisker plot for each metric. Include each group on the
% x-axis.
% INPUTS:
%   metrics_mat (double matrix): columns are samples. Rows are metrics:
%       Row 1 = length density
%       Row 2 = branch density
%       Row 3 = volume fraction
%       Row 4 = tortuosity mean
%       Row 5 = tortuosity median
%       Row 6 = tortuosity mode
%   idx (integer array): grouping index for each sample (column in
%               metrics_mat)
%   output_dir (string): output filepath
%   fstring (string): prefix for saving the figure

% Find indices in matrix corresponding to each group
g1_idx = idx==1;
g2_idx = idx==2;
g3_idx = idx==3;

% Separate each metric by grouping
g1 = metrics_mat(:,g1_idx);
g2 = metrics_mat(:,g2_idx);
g3 = metrics_mat(:,g3_idx);

% Create the titles for each plot
t_array = {'Length Density (\mum^-^2)', 'Branch Density (\mum^-^3)',...
    'Volume Fraction (a.u.)','Tortuosity Mean (a.u.)',...
    'Tortuosity Median (a.u.)','Tortuosity Mode (a.u.)'};
fig_array = {'length_density','branch_density', 'volume_fraction',...
    'tortuosity_mean','tortuosity_median','tortuosity_mode'};

% Iterate through each metric and generate bar chart
for ii = 1:size(metrics_mat,1)
    % Extract metric for each group
    x1 = g1(ii,:);
    x2 = g2(ii,:);
    x3 = g3(ii,:);
    % Create box/whisker figure
    parent=figure;
    parent.WindowState = 'maximized';
    x = [x1'; x2'; x3'];
    g = [zeros(length(x1), 1); ones(length(x2), 1); ones(length(x3), 1).*2];
    b = boxplot(x, g,'Labels',{'Group 1', 'Group 2', 'Group 3'});
    title(t_array{ii})
    set(gca, 'FontSize', 30)
    set(b,'LineWidth',4)
    % Save the output
    fig_name = append(fstring, fig_array{ii});
    fout = fullfile(output_dir, fig_name);
    saveas(gca,fout,'png');
end

end
