function [pstats] = calc_heatmap_stats(hm, regions, params, subids,...
                                        alpha, dout, fname)
%CALC_HEATMAP_STATS kolmogorov-smirnov (KS) test of distributions
% This function performs the two-sided kolmogorov-smirnov (KS) test between
% each of the experimental groups and the control groups. The null
% hypothesis of this test is that the two samples comes from the same
% distribution. The alternative hypothesis is that they come from different
% distributions.
% 
% In addition, this script measures the mean, median, and variance of each
% group for each region and parameter. These summary statistics are saved
% alongside the p-value of the KS test.
%
%   INPUTS:
%       hm (struct): contains substructures of vascular parameters
%                       metrics.[region].[parameter].[group]
%       regions (cell array): tissue regions (tiss, gyri, sulci, gm, wm)
%       params (cell array): vascular parameters (length_density,
%                            branch_density, fraction_volume, tortuosity
%                            diameter)
%       subids (cell array): subject IDs
%       alpha (double): threshold for kruskal-wallis significance
%       trend (double): threshold for a trend
%       dout (string): directory to save table
%       fname (string): filename to save table
%   OUTPUTS:
%       stats (struct): contains p-value from KS test and the summary
%           statistics for each parameter from each region:
%               - stats.[region].[parameter].p
%               - stats.[region].[parameter].
% 
%
%% Initialize workspace
% Struct for storing statistical results
pstats = struct();

% Iterate over tissue regions
for ii = 1:length(regions)
    % Iterate over parameters
    for j = 1:length(params) 
        %%% Load in the arrays from each group (ad, cte, nc)
        ad = hm.(regions{ii}).(params{j}).ad;
        cte = hm.(regions{ii}).(params{j}).cte;
        nc = hm.(regions{ii}).(params{j}).nc;

        %% Generate a linear mixed-effects model
        % The purpose is to incorporate a random effect into the analysis
        % to account for a potential correlation between outcomes on the
        % same person.

        %%% Create column vector of subject IDs for each group
        % Initialize counter for N values per group
        nval_ad = 0;
        nval_cte = 0;
        nval_nc = 0;
        % Initialize counter for number of values per subject
        nval = zeros(length(subids),1);
        % Iterate over subjects
        for s = 1:length(subids)
            % Count number of values in array
            tmp = hm.(subids{s}).(regions{ii}).(params{j});
            nval(s) = length(tmp);
            % Add N values to the respective group list
            if contains(subids(s), 'AD_')
                nval_ad = nval_ad + nval(s);
            elseif contains(subids(s), 'CTE_')
                nval_cte = nval_cte + nval(s);
            elseif contains(subids(s), 'NC_')
                nval_nc = nval_nc + nval(s);
            end
        end
        % Create column vector with subject identifier for each value
        subID_vec = zeros(sum(nval),1);
        % Start index
        s0 = 1;
        for s = 1:length(subids)
            % Update end index
            sf = s0 + nval(s) - 1;
            % Add subject ID for each value
            subID_vec(s0:sf) = ones(nval(s),1) .* s;
            % Update start index
            s0 = sf + 1;
        end
        % Separate the groups (AD, CTE, NC)
        ad_vec = subID_vec(1 : nval_ad);
        cte_vec = subID_vec(nval_ad+1 : nval_ad+nval_cte);
        nc_vec = subID_vec(nval_ad+nval_cte+1 : nval_ad+nval_cte+nval_nc);


        %%% Create a table for fitting the LME model
        % Column 1 = group label
        % Column 2 = vascular metric value
        % Create the group labels
        g_ad = repmat({'AD'},[length(ad),1]);
        g_cte = repmat({'CTE'},[length(cte),1]);
        g_nc = repmat({'HC'},[length(nc),1]);
        groups_ad_nc = vertcat(g_ad, g_nc);
        groups_cte_nc = vertcat(g_cte, g_nc);
        % Combine vascular metrics into column vector
        vmets_ad_nc = [ad; nc];
        vmets_cte_nc = [cte; nc];
        % Combine the subject IDs into column vector
        subids_ad_nc = vertcat(ad_vec, nc_vec);
        subids_cte_nc = vertcat(cte_vec, nc_vec);
        % Table (group labels, subjectID, vascular metric values)
        tbl_ad_nc = table(groups_ad_nc,subids_ad_nc,vmets_ad_nc,...
              'VariableNames',{'Groups','subID','VascularMetric'});
        tbl_cte_nc = table(groups_cte_nc,subids_cte_nc,vmets_cte_nc,...
              'VariableNames',{'Groups','subID','VascularMetric'});
               
        %%% Fit linear mixed-effects model (LME) (table)
        % Define the model:
        %   response = vascular metric array
        %   random effect (intercept/subject): subID
        %   fixed effect: Groups
        fml = 'VascularMetric ~ Groups + (1 | subID)';
        % Fit the model for AD vs. HC
        lme_ad_nc = fitlme(tbl_ad_nc,fml);
        % Fit the model for CTE vs. HC
        lme_cte_nc = fitlme(tbl_cte_nc,fml);
        
        %%% Estimates of random effects 
        [~,~,stats_ad_nc] = fixedEffects(lme_ad_nc);
        [~,~,stats_cte_nc] = fixedEffects(lme_cte_nc);
        % Store the p-value for the stats tests
        p = zeros(2,1);
        p(1) = stats_ad_nc{2,6};
        p(2) = stats_cte_nc{2,6};
        
        %%% print region/parameter/inequal. for significance/trend
        % Median values
        med_ad = median(ad);
        med_cte = median(cte);
        med_nc = median(nc);
        % Mean values
        mean_ad = mean(ad);
        mean_cte = mean(cte);
        mean_nc = mean(nc);

        %% Print the significant results
        % AD Significance
        if p(1) < alpha
            % AD Mean
            if mean_ad < mean_nc
                fprintf('SIG: mean(ad) < mean(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            elseif mean_ad == mean_nc
                fprintf('SIG: mean(ad) = mean(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            else
                fprintf('SIG: mean(ad) > mean(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            end
            % AD Median
            if med_ad < med_nc
                fprintf('SIG: median(ad) < median(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            elseif med_ad == med_nc
                fprintf('SIG: median(ad) = median(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            else
                fprintf('SIG: median(ad) > median(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            end
        end
        % CTE Significance
        if p(2) < alpha
            % AD Mean
            if mean_cte < mean_nc
                fprintf('SIG: mean(cte) < mean(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            elseif mean_cte == mean_nc
                fprintf('SIG: mean(cte) = mean(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            else
                fprintf('SIG: mean(cte) > mean(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            end
            % AD Median
            if med_cte < med_nc
                fprintf('SIG: median(cte) < median(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            elseif med_cte == med_nc
                fprintf('SIG: median(cte) = median(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            else
                fprintf('SIG: median(cte) > median(nc) for %s, %s\n',...
                        regions{ii}, params{j})
            end
        end
        
        %%% Save summary statistics
        % Save coefficients, which is the difference between the
        % experimental and control group, adjusted for potential
        % correlation between repeated measures on same subject.
        pstats.(regions{ii}).(params{j}).coef.ad_nc = stats_ad_nc{2,2};
        pstats.(regions{ii}).(params{j}).coef.cte_nc = stats_cte_nc{2,2};

        % Save the 95% confidence interval for the coefficient
        pstats.(regions{ii}).(params{j}).CI_lower.ad_nc = stats_ad_nc{2,7};
        pstats.(regions{ii}).(params{j}).CI_upper.ad_nc = stats_ad_nc{2,8};
        pstats.(regions{ii}).(params{j}).CI_lower.cte_nc = stats_cte_nc{2,7};
        pstats.(regions{ii}).(params{j}).CI_upper.cte_nc = stats_cte_nc{2,8};

        % Save the p-value for each comparison
        pstats.(regions{ii}).(params{j}).p.ad_nc = p(1);
        pstats.(regions{ii}).(params{j}).p.cte_nc = p(2);

        % Save the median, mean, & std dev for each comparison
        pstats.(regions{ii}).(params{j}).med.ad = med_ad;
        pstats.(regions{ii}).(params{j}).med.cte = med_cte;
        pstats.(regions{ii}).(params{j}).med.nc = med_nc;
        pstats.(regions{ii}).(params{j}).mean.ad = mean_ad;
        pstats.(regions{ii}).(params{j}).mean.cte = mean_cte;
        pstats.(regions{ii}).(params{j}).mean.nc = mean_nc;
        pstats.(regions{ii}).(params{j}).std.ad = std(ad);
        pstats.(regions{ii}).(params{j}).std.cte = std(cte);
        pstats.(regions{ii}).(params{j}).std.nc = std(nc);

        close all;
    end
    %% Generate a table of p-values for the region
    % Create cell array for the pairwise comparisons
    Metric = {'AD vs. HC p-value'; 'CTE vs. HC p-value';...
            'AD vs. HC coeff'; 'CTE vs. HC coeff'; 
            'AD vs. HC C.I. (lower)'; 'CTE vs. HC C.I.(lower)';...
            'AD vs. HC C.I. (upper)'; 'CTE vs. HC C.I. (upper)';...
            'AD Med.';'CTE Med.';'HC Med.';...
            'AD Std.Dev.';'CTE Std.Dev.';'HC Std.Dev.'};
    
    %{
    %%% Retrieve the p-values for each parameter
    VolumeFraction = cell2mat(struct2cell(pstats.(regions{ii}).(params{1}).p));
    LengthDensity = cell2mat(struct2cell(pstats.(regions{ii}).(params{2}).p));
    BranchDensity = cell2mat(struct2cell(pstats.(regions{ii}).(params{3}).p));
    Tortuosity = cell2mat(struct2cell(pstats.(regions{ii}).(params{4}).p));
    Diameter = cell2mat(struct2cell(pstats.(regions{ii}).(params{5}).p));
    
    %%% Retrieve coefficient and C.I. and add to metrics
    % Volume Fraction
    coef = cell2mat(struct2cell(pstats.(regions{ii}).(params{1}).coef));
    cil = cell2mat(struct2cell(pstats.(regions{ii}).(params{1}).CI_lower));
    ciu = cell2mat(struct2cell(pstats.(regions{ii}).(params{1}).CI_upper));
    VolumeFraction = [VolumeFraction; coef; cil; ciu];
    % Length Density
    coef = cell2mat(struct2cell(pstats.(regions{ii}).(params{2}).coef));
    cil = cell2mat(struct2cell(pstats.(regions{ii}).(params{2}).CI_lower));
    ciu = cell2mat(struct2cell(pstats.(regions{ii}).(params{2}).CI_upper));
    LengthDensity = [LengthDensity; coef; cil; ciu];
    % Branch Density
    coef = cell2mat(struct2cell(pstats.(regions{ii}).(params{3}).coef));
    cil = cell2mat(struct2cell(pstats.(regions{ii}).(params{3}).CI_lower));
    ciu = cell2mat(struct2cell(pstats.(regions{ii}).(params{3}).CI_upper));
    BranchDensity = [BranchDensity; coef; cil; ciu];

    %%% Retrieve the summary stats for each parameter (median, std. dev.)
    % Volume Fraction
    med = cell2mat(struct2cell(pstats.(regions{ii}).(params{1}).med));
    stdd = cell2mat(struct2cell(pstats.(regions{ii}).(params{1}).std));
    VolumeFraction = [VolumeFraction; med; stdd];
    % Length Density
    med = cell2mat(struct2cell(pstats.(regions{ii}).(params{2}).med));
    stdd = cell2mat(struct2cell(pstats.(regions{ii}).(params{2}).std));
    LengthDensity = [LengthDensity; med; stdd];
    % Branch Density
    med = cell2mat(struct2cell(pstats.(regions{ii}).(params{3}).med));
    stdd = cell2mat(struct2cell(pstats.(regions{ii}).(params{3}).std));
    BranchDensity = [BranchDensity; med; stdd];

    %%% Create the p-value table
    ptable = make_heatmap_ptable(Metric,VolumeFraction,LengthDensity,...
                                BranchDensity);
    %}
    
    % Initialize structure for storing 
    stat_struct = struct();
    % names of entries for Structure 
    param_name = {'VolumeFraction','LengthDensity','BranchDensity',...
        'Tortuosity','Diameter'};

    % Iterate over each parameter
    for j = 1:length(params)
        % Retrieve p-value
        p = cell2mat(struct2cell(pstats.(regions{ii}).(params{j}).p));
        % Retrieve rho coefficient & confidence interval
        coef = cell2mat(struct2cell(pstats.(regions{ii}).(params{j}).coef));
        cil = cell2mat(struct2cell(pstats.(regions{ii}).(params{j}).CI_lower));
        ciu = cell2mat(struct2cell(pstats.(regions{ii}).(params{j}).CI_upper));
        % Retrieve summary stats
        med = cell2mat(struct2cell(pstats.(regions{ii}).(params{j}).med));
        stdd = cell2mat(struct2cell(pstats.(regions{ii}).(params{j}).std));
        % Add summary stats to struct
        stat_struct.(param_name{j}) = [p; coef; cil; ciu; med; stdd];
    end
    
    % Create the summary statistics table
    ptable = make_heatmap_ptable(Metric,stat_struct.(param_name{1}),...
        stat_struct.(param_name{2}),stat_struct.(param_name{3}),...
        stat_struct.(param_name{4}),stat_struct.(param_name{5}));
    
    %%% Save table to CSV on a specific sheet
    % Create output filename
    table_out = fullfile(dout, fname);
    writetable(ptable, table_out, 'Sheet',regions{ii});
end

end