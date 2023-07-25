%% Test parallel batch array job on SCC
% Author: Mack Hyman
% Detailed Description
%{
This script performs the following:
- pass in argument from bash script
- print argument to command window
%}
clear; clc; close all;

%% Print the batch_idx
batch_idx = getenv('$SGE_TASK_ID');
sprintf(num2str(batch_idx))

%% Input variable from bash script
%{
%%% Define subject ID based upon the subid_idx from the bash script
subid_lst = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424',...
             'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_8572', 'CTE_7126',...
             'NC_21499', 'NC_6047', 'NC_6839', 'NC_6974', 'NC_7597',...
             'NC_8095', 'NC_8653'};
% subid_idx is an input parameter
subid = subid_lst{subid_idx};

%%% Std. Dev. for gaussian filter
% Define Gaussian sigma array based upon the gauss_idx from the bash script
% The standard sigma array = [1, 3, 5]
% Medium vessel sigma array = [7, 9, 11]
% Large vessel sigma array = [13, 15, 17]
sigma_lst = [1, 3, 5;...
			7, 9, 11;...
			13, 15, 17];
gsigma = sigma_lst(gauss_idx,:);

%%% Print output
sprintf(subid)
sprintf(gsigma)
%}