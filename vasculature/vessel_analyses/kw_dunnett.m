function [p] = kw_dunnett(ad, cte, nc)
%KW_DUNNETT Kruskal-Wallis followed by multi comparisons (Dunnett's)
%   This is an example function of how to compare multiple experimental
%   groups to the same control group. The Kruskal-Wallis tests whether each
%   dataset comes from the same distribution. The Dunnett's multiple
%   comparisons method corrects for multiple comparisons between multiple
%   experimental groups with the same control group.
% INPUTS:
%   ad (double column vector [n_ad,1]): the values from the AD group
%   cte (double column vector [n_ad,1]): the values from the CTE group
%   nc (double column vector [n_ad,1]): the values from the NC group
% OUTPUTS
%   p (double array): the p-values for the two experimental groups compared
%       to the control group. p(1) is the p-value for AD. p(2) is the
%       p-value for CTE.

%% Initialize variables
% array for storing P-value
p = zeros(2,1);

%%% Create labels for each vector of [experimental; control]
[nc_group{1:size(nc)}] = deal('nc');
[ad_group{1:size(ad)}] = deal('ad');
[cte_group{1:size(cte)}] = deal('cte');
% Combine all group labels into a single array. This is used in teh
% kruskal-wallis test for assigning each data point to a group
group = [ad_group, cte_group, nc_group]';

%%% Concatenate the experimental & control vectors into a matrix
% This is the input data for the kruskal-wallis test
ad_cte_nc = [ad; cte; nc];

%% Call the Kruskal Wallis Test
[~,~,stats] = kruskalwallis(ad_cte_nc, group);
% Multiple Comparisons w/ Dunnett correction
[results,~,~,~] = multcompare(stats,...
    'ControlGroup',3,...
    'CriticalValueType','dunnett');
% Extract the p-values for each comparison
p(1) = results(1,6);
p(2) = results(2,6);


end