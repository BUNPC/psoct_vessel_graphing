function yaxis_scatter(c, var, n)
%YAXIS_SCATTER vertical scatter plot for values within the same group
%   add a horizontal line for the mean value of each group.
% INPUTS:
%   c (double array [n,3]) - color of dots and lines
%   var (struct) - data for y-axis scatter plot. Each field in the struct
%       contains the data for a particular x-axis offset plot
%   n (double) - x-axis position

% Example
%{
data = readtable('PMI_stat_analysis.xlsx');
Category = ['AD', 'CTE', 'NC'];

N_def = data.N_def;
cat = data.Case;
PMI = data.PMI;
stage = data.Stage;

var1 = N_def(find(PMI<24&cat==1)); %AD
var2 = N_def(find(PMI<24&cat==2)); %CTE
var3 = N_def(find(PMI<24&cat==3)); % NC
var4 = N_def(find(PMI<24&cat==3&stage==4)); % NC Stage 4

figure(1); hold on;
c = [0.4940 0.1840 0.5560];
subplot_scat_bar(c, var1, 1)
c = [0.4660 0.6740 0.1880];
subplot_scat_bar(c, var2, 2)
c = [0.8500 0.3250 0.0980];
subplot_scat_bar(c, var3, 3)
c = [108/255 46/255 18/255];

function subplot_scat_bar
scatter(3*ones(length(var4)), var4, 'MarkerFaceColor',c, 'MarkerEdgeColor',c)
hold off
set(gca, 'FontSize', 18, 'FontWeight', 'bold');
xticklabels({'', 'AD', 'CTE', 'NC'});
ylabel('N_{def}/mm^{2}')
end
%}

%% Initialize figure
hdl = figure;
hdl.WindowState = 'maximized';
hold on;

%% Iterate through the groups in the struct

% Retrieve field names
f =  fields(var);


% Scatter plot of values within group
scatter(n*ones(length(var)), var, 'MarkerFaceColor',c, 'MarkerEdgeColor',c)

% Add horizontal line for mean value
mean_var = mean(var);
plot([n-0.4 n+0.4], [mean_var, mean_var], 'LineWidth', 3, 'Color',c);
end