%% Main script for kmeans, PCA, and graphical data analyses on CAA cases
% This was initially written by James Edelist for the Spring 2024 UROP.
% It has since been modified by Mack Hyman.
clear; clc; close all;

%% Managing environment and paths (James' laptop)
%{
% Add path to wokring directory and .mat data struct
addpath('/Users/jamieedelists/Desktop/spring_UROP_analyses');
addpath('/Users/jamieedelists/Desktop/spring_UROP_analyses/violin/');
% Load the data struct
struct_subs = load(['/Users/jamieedelists/Desktop/spring_UROP_analyses' ...
    '/caa_crop300_rmloops_metrics.mat']);
%}

%% Add top-level directory
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
output_dir = '/projectnb/npbssmic/ns/CAA/metrics/caa_entire_volume_loops';

% Load the data struct
caa_path = '/projectnb/npbssmic/ns/CAA/';
metrics = load('/projectnb/npbssmic/ns/CAA/metrics/caa_loops_metrics.mat');
metrics = metrics.metrics;

%% Sort metrics struct into array for each metric
% Retrieve subject IDs from metrics struct
subjects = fieldnames(metrics);
% Initialize/preallocate arrays to store metrics for each subject
len_den = zeros(1, numel(subjects));
len_den_skel = zeros(2,numel(subjects));
branch_den = zeros(1, numel(subjects));
frac_vol = zeros(1, numel(subjects));
tort_mean = zeros(1, numel(subjects));
tort_med = zeros(1, numel(subjects));
tort_mode = zeros(1, numel(subjects));
% Placing metrics into seperate arrays for concatenation (except tort_array)
% Extract metrics for each subject
for ii = 1:numel(subjects)
    len_den(ii) = metrics.(subjects{ii}).length_density;
    branch_den(ii) = metrics.(subjects{ii}).branch_density;
    frac_vol(ii) = metrics.(subjects{ii}).fraction_volume;
    tort_mean(ii) = metrics.(subjects{ii}).tortuosity_mean;
    tort_med(ii) = metrics.(subjects{ii}).tortuosity_med;
    tort_mode(ii) = metrics.(subjects{ii}).tortuosity_mode;
    % Load skeleton length density
    len_den_skel(:,ii) = metrics.(subjects{ii}).length_density_skeleton;
end
% Combine metrics into a matrix. Rows = metrics. Columns = samples.
metrics_mat = vertcat(len_den, branch_den, frac_vol, tort_mean,...
                    tort_med, tort_mode);

%% Compare length densities between skeleton and graph

%%% Plot the upper bounds of skeleton
parent = figure;
parent.WindowState = 'maximized';
% Retrieve upper bounds
ld_skel_upper = len_den_skel(2,:);
% Set the width of the overlaid bars
w1 = 0.9;
bar(1:size(subjects,1), ld_skel_upper, w1, 'FaceColor',[0.2 0.2 0.5]);

%%% Overlay the length density from the graph
hold on;
% Set the width of the overlaid bars
w2 = 0.6;
bar(1:size(subjects,1), len_den, w2,'FaceColor',[1 0 0]);

%%% Overlay the lower bounds of skeleton
% Set the width of the overlaid bars
w3 = 0.3;
ld_skel_lower = len_den_skel(1,:);
bar(1:size(subjects,1), ld_skel_lower,w3,'FaceColor',[0 0.7 0.7]);

% Labels and titles
xticklabels(subjects);
set(gca,'XTickLabelRotation',30);
set(gca,'TickLabelInterpreter', 'none');
set(gca, "FontSize", 30);
legend({'Skeleton Upper Bound', 'Graph', 'Skeleton Lower Bound'});
title('Length Density Comparisons')
ylabel('Total Length (\mum) / Volume (\mum^3)')


%%% Save the chart
fig_name = 'length_density_comparisons';
fout = fullfile(output_dir, fig_name);
saveas(gcf,fout,'png');

%% Create bar charts (x-axis location for each subject)
% Cell array for name of each metric
metrics_names = {'Length Density', 'Branch Density', 'Fraction Volume',...
    'Tortuosity Mean', 'Tortuosity Median', 'Tortuosity Mode'};
fig_array = {'length_density','branch_density', 'volume_fraction',...
    'tortuosity_mean','tortuosity_median','tortuosity_mode'};
ylabels = {'Length (\mum) / \mum^3','Branches / \mum^3','(a.u.)','(a.u.)','(a.u.)','(a.u.)'};
bar_chart(metrics_mat, metrics_names, ylabels,...
    fig_array, subjects, output_dir);

%% k-means clustering on all metrics (except tortuosity)
%%% k-means
% Number of clusters (three disease stages)
k = 3;
% kMeans matrix variables (row = case, column = variable)
kmeans_x = [len_den; branch_den; frac_vol; tort_mean; tort_med; tort_mode]';
% Perform k-means clustering with data
[idx, c, sumd] = kmeans(kmeans_x, k);

%%% Create box/whisker plots kmeans1 (without PCA)
% Create box/whisker plot
fstring = 'kmeans1_';
group_box_whisker(metrics_mat, idx, output_dir, fstring)

%% Principal components analysis

%%% Perform PCA on the data
[coeff, score, latent, tsquared, explained] = pca(kmeans_x);

%%% Choose the number of components to keep based on explained variance
% Find components that explain at least 85% of the variance
explained_sum = cumsum(explained);
% Find components that explain at least (threshold)% of the variance
threshold = 95;
numComponents = find(explained_sum >= threshold, 1, 'first');
% Reduce data dimensions based on the number of principal componets (first
% row or first 2 rows)
x_reduced = score(:, 1:numComponents);

%% Perform k-means clustering on reduced data (after PCA)
close all;
%%% k-means with PCA data
[idx2, C2] = kmeans(x_reduced, k);

%%% Graphing the kMeans clusters of principal components
% WITHOUT COLORING THE BACK OF THE GRAPH INTO PARTITIONS OF CLUSTERS
% HOWEVER THIS WORKS IF THERE ARE 1, 2, or 3 PRINCIPLE COMPONENTS
% Check the number of principal components and visualize accordingly
if numComponents == 1
    figure;
    scatter(x_reduced(:,1), ones(size(x_reduced, 1), 1), 10, idx, 'filled');
    xlabel('Principal Component 1');
    title('PCA followed by k-Means Clustering (1D)');
elseif numComponents == 2
    figure;
    gscatter(x_reduced(:,1), x_reduced(:,2), idx);
    xlabel('Principal Component 1');
    ylabel('Principal Component 2');
    title('PCA followed by k-Means Clustering (2D)');
elseif numComponents == 3
    figure;
    scatter3(x_reduced(:,1), x_reduced(:,2), x_reduced(:,3), 10, idx, 'filled');
    xlabel('Principal Component 1');
    ylabel('Principal Component 2');
    zlabel('Principal Component 3');
    view(3);
    title('PCA followed by k-Means Clustering (3D)');
end

%%% Generate barcharts with k-means (after PCA)
fstring = 'kmeans2_';
group_box_whisker(metrics_mat, idx2, output_dir, fstring)

%%% Create and print table of cluster assignments
t = table(subjects, idx, idx2);
disp(t)

%% GRAPH THE KMEANS CLUSTERS OF TWO PRINCIPAL COMPONENTS AND COLORS BACKGROUND OF PLOT
% INTO PARTITIONS BASED ON CLUSTERS OF THE 2 PRINCIPAL COMPONENTS
% ONLY WORKS IF THERE ARE TWO PRINTIPAL COMPONENTS (NOT GENERAL)
% https://www.mathworks.com/help/stats/kmeans.html

% Initialize full size window
parent=figure;
parent.WindowState = 'maximized';

% Create linear space for meshgrid of the kmeans
x_range = linspace(min(x_reduced(:,1)), max(x_reduced(:,1)), 500);
y_range = linspace(min(x_reduced(:,2)), max(x_reduced(:,2)), 500);
[X1, X2] = meshgrid(x_range, y_range);
grid_X = [X1(:), X2(:)];
[grid_clusters, ~] = kmeans(grid_X, k, 'MaxIter', 1, 'Start', C2);

% Plot the grid with cluster colors
gscatter(grid_X(:,1), grid_X(:,2), grid_clusters,...
    [0,0.75,0.75;0.75,0,0.75;0.75,0.75,0], '..');
hold on;

% Overlay grid with PCA data points
plot(x_reduced(:,1), x_reduced(:,2), 'k*', ...
    'MarkerSize', 15, 'LineWidth', 10);
title('PCA followed by K-Means Clustering (2D)');
xlabel('Principal Component 1');
ylabel('Principal Component 2');
[~,lgd] = legend('Group 2','Group 1','Group 3','Data','Location','SouthEast',...
    'FontSize',15);
% Set marker size in legend
ob = findobj(lgd, 'type', 'line');
set(ob, 'Markersize',20);
set(gca, 'FontSize', 30)
hold off;

% Save output
file_out = fullfile(output_dir, 'kmeans_after_pca');
saveas(gcf, file_out, 'png');

%% Box and Whiskers plot for length_den, branch_den, frac_vol

%%% Set font size and line width for all
font_size = 25;
line_width = 4;

%%% length_density
parent=figure;
parent.WindowState = 'maximized';
% Draw Boxplot
boxplot(len_den', 'Labels', {'All CAA Subjects'});
title('Length Density');
ylabel('Length (\mum) / \mum^3');
set(gca, 'FontSize', font_size);   %increase font size to 30
set(findall(gca, 'Type', 'line'), 'LineWidth', line_width); %increase line width
% Save figure
file_out = fullfile(output_dir, 'box_whisker_length_density_all_except_caa22occ');
saveas(gcf, file_out, 'png');

%%% Subplot for branch_density
parent=figure;
parent.WindowState = 'maximized';
% Draw Boxplot
boxplot(branch_den', 'Labels', {'All CAA Subjects'});
title('Branch Density');
ylabel('Branches / \mum^3');
set(gca, 'FontSize', font_size);
set(findall(gca, 'Type', 'line'), 'LineWidth', line_width);
% Save figure
file_out = fullfile(output_dir, 'box_whisker_branch_density_all_except_caa22occ');
saveas(gcf, file_out, 'png');

%%% Subplot for fraction_volume
parent=figure;
parent.WindowState = 'maximized';
% Draw Boxplot
boxplot(frac_vol', 'Labels', {'All CAA Subjects'});
title('Fraction Volume');
ylabel('(a.u.)');
set(gca, 'FontSize', font_size);
set(findall(gca, 'Type', 'line'), 'LineWidth', line_width);
% Save figure
file_out = fullfile(output_dir, 'box_whisker_fraction_volume_all_except_caa22occ');
saveas(gcf, file_out, 'png');

% Main title for entire figure
% sgtitle('All Samples Combined');

%% Violin plots for tortuosity arrays
%{
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
set(gca, 'FontSize', 24); % increase font size
%}

%% Function to create bar chart for each metric (except tortuosity)
function bar_chart(metrics_mat, metrics_names, ylabels, fig_array,...
    subjects, output_dir)
% BAR_CHART bar chart to visualize value for each subject
% INPUTS:
%   metrics_mat (double matrix): columns are samples. Rows are metrics:
%       Row 1 = length density
%       Row 2 = branch density
%       Row 3 = volume fraction
%   metrics_names (cell array): name of each metric
%   fig_array (cell array): names of each figure to save
%   ylabels (cell array): y-axis labels/units
%   subjects (cell array): cell of subject ID
%   output_dir (string): output filepath
%   fstring (string): prefix for saving the figure

%%% Iterate through each metric (rows in metrics_mat)
for ii = 1:size(metrics_mat,1)
    % Extract y-axis values for this metric
    y = metrics_mat(ii,:);
    % Create a bar-chart
    parent=figure;
    parent.WindowState = 'maximized';
    bar(1:size(subjects,1), y);
    % Labels and titles
    xticklabels(subjects);
    ylabel(ylabels{ii});
    title(metrics_names(ii));
    set(gca,'XTickLabelRotation',30);
    set(gca,'TickLabelInterpreter', 'none');
    set(gca, "FontSize", 30);
    % Save the chart
    fig_name = string(append(fig_array(ii), '_barchart'));
    fout = fullfile(output_dir, fig_name);
    saveas(gcf,fout,'png');
end

end
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