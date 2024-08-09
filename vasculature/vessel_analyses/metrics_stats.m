function [pstats] = metrics_stats(metrics, regions, params,...
                        alpha, trend, dout, fname)
%METRICS_ANOVA perform Kruskal-Wallis on the subset of metrics
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
%   correction.
%   INPUTS:
%       metrics (struct): contains substructures of vascular parameters
%                       metrics.[region].[parameter].[group]
%       regions (cell array): tissue regions (tiss, gyri, sulci, gm, wm)
%       params (cell array): vascular parameters (length_density,
%                              branch_density, fraction_volume, tortuosity)
%       alpha (double): threshold for kruskal-wallis significance
%       trend (double): threshold for a trend
%       ntests (uint): number of comparisons. this is used for multiple
%                       comparisons corrections
%       dout (string): directory to save table
%       fname (string): filename to save table
%   OUTPUTS:
%       stats (struct): contains p-value and test decision for either ANOVA
%           or kurskal-wallis test for each parameter from each region.
%               p-value: stats.[region].[parameter].p
%
% TODO:
%   - change kstest2 to the LME model

%% Initialize workspace
% Struct for storing statistical results
pstats = struct();

%% Iterate over tissue regions
for ii = 1:length(regions)
    %%% Iterate over parameters
    for j = 1:length(params) 
        %% Concatenate groups into matrix for pair-wise comparison
        % Load in the arrays from each group (ad, cte, nc)
        ad = metrics.(regions{ii}).(params{j}).ad;
        cte = metrics.(regions{ii}).(params{j}).cte;
        nc = metrics.(regions{ii}).(params{j}).nc;
        
        %%% multiple comparisons (kruskal-wallis)
        p = zeros(2,1);
        % Create labels for each vector of [experimental; control]
        [nc_group{1:size(nc)}] = deal('nc');
        [ad_group{1:size(ad)}] = deal('ad');
        [cte_group{1:size(cte)}] = deal('cte');
        group = [ad_group, cte_group, nc_group]';
        % Concatenate the experimental & control into vector
        ad_cte_nc = [ad; cte; nc];
        % Call the Kruskal Wallis Test
        [~,~,stats] = kruskalwallis(ad_cte_nc, group);
        % Multiple Comparisons w/ Dunnett correction
        [results,~,~,~] = multcompare(stats,...
            'ControlGroup',3,...
            'CriticalValueType','dunnett');
        % Extract the p-values for each comparison
        p(1) = results(1,6);
        p(2) = results(2,6);

        
        %%% experimental/healthy comparisons for statistic
        ad_mean_tf = mean(ad) < mean(nc);
        ad_med_tf = median(ad) < median(nc);
        cte_mean_tf = mean(cte) < mean(nc);
        cte_med_tf = median(cte) < median(nc);

        %%% print region/parameter/inequal. for significance/trend
        % AD Significance
        if p(1) < alpha
            % AD Mean
            if ad_mean_tf == 1
                fprintf('SIG: mean(ad) < mean(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            else
                fprintf('SIG: mean(ad) > mean(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            end
            % AD Median
            if ad_med_tf == 1
                fprintf('SIG: median(ad) < median(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            else
                fprintf('SIG: median(ad) > median(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            end
        end
        % AD Trend
        if (alpha < p(1)) && (p(1) < trend)
            % AD Mean
            if ad_mean_tf == 1
                fprintf('TREND: mean(ad) < mean(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            else
                fprintf('TREND: mean(ad) > mean(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            end
            % AD Median
            if ad_med_tf == 1
                fprintf('TREND: median(ad) < median(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            elseif ad_med_tf == 0
                fprintf('TREND: median(ad) > median(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            end
        end
        % CTE Significance
        if p(2) < alpha
            % AD Mean
            if cte_mean_tf == 1
                fprintf('SIG: mean(cte) < mean(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            else
                fprintf('SIG: mean(cte) > mean(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            end
            % AD Median
            if ad_med_tf == 1
                fprintf('SIG: median(cte) < median(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            else
                fprintf('SIG: median(cte) > median(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            end
        end
        % CTE Trend
        if (alpha < p(2)) && (p(2) < trend)
            % CTE Mean
            if cte_mean_tf == 1
                fprintf('TREND: mean(cte) < mean(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            else
                fprintf('TREND: mean(cte) > mean(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            end
            % CTE Median
            if cte_med_tf == 1
                fprintf('TREND: median(cte) < median(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            else
                fprintf('TREND: median(cte) > median(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            end
        end


        % Save the p-value for each comparison
        pstats.(regions{ii}).(params{j}).p.ad_nc = p(1);
        pstats.(regions{ii}).(params{j}).p.cte_nc = p(2);
        close all;
    end

    %% Generate a table of p-values for the region
    % Create cell array for the pairwise comparisons
    Pairs = {'AD vs. HC'; 'CTE vs. HC'};
    % Retrieve the p-values for each parameter
    LengthDensity = cell2mat(struct2cell(pstats.(regions{ii}).(params{1}).p));
    BranchDensity = cell2mat(struct2cell(pstats.(regions{ii}).(params{2}).p));
    VolumeFraction = cell2mat(struct2cell(pstats.(regions{ii}).(params{3}).p));
    TortOutliers = cell2mat(struct2cell(pstats.(regions{ii}).(params{4}).p));
    Tort = cell2mat(struct2cell(pstats.(regions{ii}).(params{5}).p));
    Diameter = cell2mat(struct2cell(pstats.(regions{ii}).(params{6}).p));
    % Create the p-value table
    ptable = make_ptable(Pairs, LengthDensity, BranchDensity,...
                VolumeFraction, TortOutliers,Tort,Diameter);
    
    %%% Save table to CSV on a specific sheet
    % Create output filename
    table_out = fullfile(dout, fname);
    writetable(ptable, table_out, 'Sheet',regions{ii});
end
end