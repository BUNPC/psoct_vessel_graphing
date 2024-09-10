function [p, coef, lower, upper] =...
    lme_stats(exp_values, exp_sub_id, cntrl_values, cntrl_sub_id)
%LME_STATS Linear Mixed Effects (LME) model
% This function creates an LME model and then tests the fixed effects
% (diseased vs. control) while accounting for random effects within groups.
% This method accounts for multiple observations (samples) from the same
% subject.
%
%   INPUTS:
%       exp_values (vector): concatenated vector of all values for all
%           subjects in the experimental group
%       exp_sub_id (vector): subject identifier for each value in vector
%           "exp_values". This is to identify which values belong to to
%           which subject. For example [1, 1, 1, ..., 2,2,2...,N,N,N] where
%           N represents the total number of subjects in the group.
%       cntrl_values (vector): same as "exp_values" except for the control
%           group.
%       cntrl_sub_id (vector): same as "exp_sub_id" except for the control
%           group.
%   OUTPUTS:
%       p (float): p-value from the fixed effects test of the LME model
% 

%% Generate a linear mixed-effects model
% The purpose is to incorporate a random effect into the analysis
% to account for a potential correlation between outcomes on the
% same person.

%%% Create a table for fitting the LME model
% Column 1 = group label ('experimental' or 'control'
% Column 2 = value of each observation
% Create the group labels
g_exp = repmat({'EXP'},[length(exp_values),1]);
g_cntrl = repmat({'CNTRL'},[length(cntrl_values),1]);
grp_lbl = vertcat(g_exp, g_cntrl);
% Combine observations into column vector
obs = [exp_values; cntrl_values];
% Combine the subject IDs into column vector
subids = vertcat(exp_sub_id, cntrl_sub_id);
% Table (group labels, subjectID, vascular metric values)
tbl_ad_nc = table(grp_lbl,subids,obs,...
                  'VariableNames',{'Groups','subID','Observation'});
       
%%% Fit linear mixed-effects model (LME) (table)
% Define the model:
%   response = measurement (Observation)
%   random effect (intercept/subject): subject identifier (subID)
%   fixed effect: group (experimental vs. control) (Groups)
fml = 'Observation ~ Groups + (1 | subID)';
% Fit the model for AD vs. HC
lme_model = fitlme(tbl_ad_nc,fml);

%% Estimates of fixed effects (statistical test) of LME model
[~,~,stats] = fixedEffects(lme_model);
% Estimated coefficient value (which is the difference between the
% experimental and control group, adjusted for potential correlation
% between repeated measures on same subject). The coefficient is the
% difference between the mean of group 2 (control) minus the mean of group
% 1 (experimental).
coef = stats{2,2};
% Lower limit of 95% conf. interval for fixed effect coefficient
lower = stats{2,7};
% Upper limit of 95% conf. interval for fixed effect coefficient
upper = stats{2,8};
% Retrieve the p-value for the stats tests
p = stats{2,6};
end