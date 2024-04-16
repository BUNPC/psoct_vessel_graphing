function bar_chart(ad_param, cte_param, nc_param, subjects, mpath,...
    tstring, xstring, ystring, fout)
%BAR_CHART Create bar chart for each metric
% Sort the input parameters into a categorical variable, generate bar
% chart, and save the bar chart.
% INPUTS:
%   ad_param (array): array of values for the parameter in the AD group
%   cte_param (array): array of values for the parameter in the CTE group
%   nc_param (array): array of values for the parameter in the NC group
%   subjects (cell array): list of cells of the subject IDs
%   mpath (string): output directory for the barchart
%   tstring (string): title of bar chart
%   xstring (string): x-axis label of bar chart
%   ystring (string): y-axis label of bar chart
%   fout (string): output filename for barchart

%%% Sort the input parameters
parameter = vertcat(ad_param, cte_param, nc_param);

%%% Initialize the bar chart
fh = figure();
fh.WindowState = 'maximized';
% Create x-axis labels according subject ID list (AD, CTE, NC)
x = categorical(subjects);
% Create bar chart
b = bar(x, parameter);
title(tstring)
xlabel(xstring)
ylabel(ystring)
set(gca, 'FontSize', 30)
set(gca, 'TickLabelInterpreter','none')

%%% Set color of bars for each group (AD, CTE, NC)
% Find number of each group
n_ad = length(find(contains(subjects,'AD')));
n_nc = length(find(contains(subjects,'NC')));
% Set face color so colors can be adjusted numerically
b.FaceColor = 'flat';
% Set the colors for AD (purple)
b.CData(1:n_ad,:) = repmat([0.5, 0, 0.5],n_ad,1);
% Set the colors for NC (green)
b.CData(end-(n_nc-1):end, :) = repmat([0,0.5,0],n_nc,1);

%%% Save output
mout = fullfile(mpath, fout);
saveas(gca, mout, 'png')

end