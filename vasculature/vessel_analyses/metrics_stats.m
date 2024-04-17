function [pstats] = metrics_stats(metrics, regions, params,...
                        th_trend, th_sigdif)
%METRICS_ANOVA perform ANOVA on the subset of metrics
%   Parse the "metrics" struct, separate each into its constituent vascular
%   metrics, test for normality (lillietest). Then, perform either ANOVA
%   (for the normally distributed parameters) or the kurskal-wallis test
%   for the non-normally distributed parameters.
%   INPUTS:
%       metrics (struct): contains substructures of vascular parameters
%                       metrics.[region].[parameter].[group]
%       regions (cell array): tissue regions (tiss, gyri, sulci, gm, wm)
%       params (cell array): vascular parameters (length_density,
%                              branch_density, fraction_volume, tortuosity)
%       groups (cell array): groups to compare (ad, cte, nc)
%       th_trend (double): threshold for a trend
%       th_sigdif (double): threshold for a significant difference
%   OUTPUTS:
%       stats (struct): contains p-value and test decision for either ANOVA
%           or kurskal-wallis test for each parameter from each region.
%               p-value: stats.[region].[parameter].p
%               test decision: stats.[region].[parameter].h

%% Test each group for normality
pstats = struct();
% Iterate over tissue regions
for ii = 1:length(regions)
    % Iterate over parameters
    for j = 1:length(params)
        %% Skip tortuosity for the ratio metrics
        if strcmp(params{j},'tortuosity') && contains(regions{ii},'sulci_')
            continue
        else
            %% Check for normality
            % Load in the arrays from each group (ad, cte, nc)
            ad = metrics.(regions{ii}).(params{j}).ad;
            cte = metrics.(regions{ii}).(params{j}).cte;
            nc = metrics.(regions{ii}).(params{j}).nc;
            % Check for normality (lillietest)
            lt = zeros(1,3);
            lt(1) = lillietest(ad);
            lt(2) = lillietest(cte);
            lt(3) = lillietest(nc);
    
            %% Concatenate groups into matrix for unbalanced comparison
            % Create labels arrays
            ad_label = repmat("AD",1,length(ad));
            cte_label = repmat("CTE",1,length(cte));
            nc_label = repmat("NC",1,length(nc));
            label_array = [ad_label, cte_label, nc_label];
    
            % Combine each group into single 1D array
            metric_array = [ad', cte', nc'];
            
            % ANOVA: all distributions are normal
            if ~all(lt)
                % Save result of normality
                pstats.(regions{ii}).(params{j}).normality = true;
                % Perform ANOVA
                [p, ~, stats] = anova1(metric_array, label_array);
            % kruskal-wallis test: at least one distribution not normal
            else
                % Save result of normality
                pstats.(regions{ii}).(params{j}).normality = false;
                % Perform kruskal-wallis test
                [p, ~, stats] = kruskalwallis(metric_array, label_array);
            end
            % Examine test statistic of each group comparison
            results = multcompare(stats);
            % Save the p-value for each comparison
            pstats.(regions{ii}).(params{j}).p.group = p;
            pstats.(regions{ii}).(params{j}).p.ad_cte = results(1,end);
            pstats.(regions{ii}).(params{j}).p.ad_nc = results(2,end);
            pstats.(regions{ii}).(params{j}).p.cte_nc = results(3,end);
            % Print the region/parameter for p-values below threshold
            if p <= th_sigdif
                fprintf('Significant difference within the %s for %s\n',...
                    regions{ii}, params{j})
            elseif p <= th_trend
                fprintf('Trend within the %s for %s\n',...
                    regions{ii}, params{j})
            end
            close all;
        end
    end
end

end