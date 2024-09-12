%% Test the function "lme_stats" using data from Anna
% This function was writting for testing Anna's experimental data. This
% test script was written to test the function.
clear; clc;


%% Import the spreadsheet of data

% Import spreadsheet of mus values
data = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/' ...
        'mus_ratio_ret_age_Sep_2024.xls'];
T = readtable(data);

% Separate table into subID, group, mus
subid = table2array(T(:,1));
grp = table2array(T(:,3));
mus = table2array(T(:,4));

% Extracat row identifier for each group by "Category" entry in table
grp1 = find(grp==1);
grp2 = find(grp==2);
grp3 = find(grp==3);

%%% Extract subID, mus for each group
% Extract subject ID for each group
grp1_subid = subid(grp1,:);
grp2_subid = subid(grp2,:);
grp3_subid = subid(grp3,:);

% Extract subject ID for each group
grp1_mus = mus(grp1,:);
grp2_mus = mus(grp2,:);
grp3_mus = mus(grp3,:);

%% Call lme_stats function on a subset of data
% Compare group 2 and 3
[p, coef, lower, upper] =...
    lme_stats(grp1_mus,grp1_subid,grp2_mus,grp2_subid);