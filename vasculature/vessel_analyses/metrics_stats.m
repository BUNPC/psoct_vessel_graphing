function [pstats] = metrics_stats(metrics, regions, params, groups,...
                        alpha, trend)
%METRICS_ANOVA perform ANOVA on the subset of metrics
%   Parse the "metrics" struct, separate each into its constituent vascular
%   metrics. Given the small sample size, it is unlikely that any of these
%   samples follow a normal distribution. Therefore, the kruskal-wallis
%   test will be used to test the null hypothesis for each metric, that all
%   samples originate from the same distribution. This initial test will
%   only indicate whether at least one sample is from a different
%   distribution, but it does not report which.
% 
%   For the metrics that reject the null hypothesis of the KW test, this
%   function will then perform a non-parametric multiple comparisons test.
%   In this case, it will be the Wilcoxon rank sum test with bonferoni
%   correction. The number of hypotheses is calculated from the number of
%   groups (factorial(N-1)).
%   INPUTS:
%       metrics (struct): contains substructures of vascular parameters
%                       metrics.[region].[parameter].[group]
%       regions (cell array): tissue regions (tiss, gyri, sulci, gm, wm)
%       params (cell array): vascular parameters (length_density,
%                              branch_density, fraction_volume, tortuosity)
%       groups (cell array): groups to compare (ad, cte, nc)
%       alpha (double): threshold for kruskal-wallis significance
%       trend (double): threshold for a trend
%   OUTPUTS:
%       stats (struct): contains p-value and test decision for either ANOVA
%           or kurskal-wallis test for each parameter from each region.
%               p-value: stats.[region].[parameter].p
%               test decision: stats.[region].[parameter].h
%
%% Initialize workspace
% Struct for storing statistical results
pstats = struct();
%%% Calculate the bonferroni-corrected significance level
n = size(groups,2);
syms k
% Find the sum of the series k from k = 1 to (N groups - 1)
m = double(symsum(k,k,1,n-1));

% bonferroni correction
alpha = alpha ./ m;


%% Iterate over tissue regions
for ii = 1:length(regions)
    % Iterate over parameters
    for j = 1:length(params)
        %% Skip tortuosity for the ratio metrics
        if strcmp(params{j},'tortuosity') && contains(regions{ii},'sulci_')
            continue
        else    
            %% Concatenate groups into matrix for unbalanced comparison
            % Load in the arrays from each group (ad, cte, nc)
            ad = metrics.(regions{ii}).(params{j}).ad;
            cte = metrics.(regions{ii}).(params{j}).cte;
            nc = metrics.(regions{ii}).(params{j}).nc;
            % Create labels arrays
            ad_label = repmat("AD",1,length(ad));
            cte_label = repmat("CTE",1,length(cte));
            nc_label = repmat("NC",1,length(nc));
            label_array = [ad_label, cte_label, nc_label];
    
            % Combine each group into single 1D array
            metric_array = [ad', cte', nc'];
            
            %%% kruskal-wallis test: at least one distribution not normal
            % Perform kruskal-wallis test
            [p, ~, ~] = kruskalwallis(metric_array, label_array);
            % Save the result
            pstats.(regions{ii}).(params{j}).p.group = p;
            
            %%% pair-wise comparisons
            p = ones(3,1);
            % Use kolmogorov-smirnov test for tortuosity distributions
            if strcmp(params{j},'tortuosity')
                [~,p(1)] = kstest2(ad, cte,'Alpha',alpha);
                [~,p(2)] = ranksum(ad, nc,'Alpha',alpha);
                [~,p(3)] = ranksum(cte, nc,'Alpha',alpha);
            % Use Wilcoxon rank sum test for all others
            else
                p(1) = ranksum(ad, cte);
                p(2) = ranksum(ad, nc);
                p(3) = ranksum(cte, nc);
            end
            % Print region/parameter for p-values below alpha/threshold
            if any(p < alpha)
                fprintf('Significant difference within the %s for %s\n',...
                        regions{ii}, params{j})
            elseif any(p < trend)
                fprintf('Trend within the %s for %s\n',...
                        regions{ii}, params{j})
            end

            % Save the p-value for each comparison
            pstats.(regions{ii}).(params{j}).p.ad_cte = p(1);
            pstats.(regions{ii}).(params{j}).p.ad_nc = p(2);
            pstats.(regions{ii}).(params{j}).p.cte_nc = p(3);
            close all;
        end
    end
end

end