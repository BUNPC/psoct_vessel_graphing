%% Test parallel batch array job on SCC
% Author: Mack Hyman
% Detailed Description
%{
This script performs the following:
- pass in argument from bash script
- print argument to command window
%}
clear; clc; close all;

%% Test first part of main script

% Path to top-level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_10T/';
% Subfolder containing data
subdir = '/dist_corrected/volume/';
% Filename to parse (this will be the same for each subject)
fname = 'ref_4ds_norm_inv';
% filename extension
ext = '.tif';

% Complete subject ID list for Ann_Mckee_samples_10T
subid = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_8572', 'CTE_7126',...
         'NC_21499', 'NC_6047', 'NC_6839', 'NC_6974', 'NC_7597',...
         'NC_8095', 'NC_8653'};

% Gaussian sigma arrays:
% Small vessel sigma array = [1, 3, 5]
% Medium vessel sigma array = [7, 9, 11]
% Large vessel sigma array = [13, 15, 17]
sigmas = [1, 3, 5; 7, 9, 11; 13, 15, 17];

%%% Create cell array of subject ID and sigma for job array on the SCC 
nrow = length(subid)*size(sigmas,2);
nsigma = size(sigmas,2);
sub_sigma = cell(length(subid).*size(sigmas,2), 2);
idx = 1;
% Fill sub_sigma cell array with each sigma array for each subject
for i = 1:3:nrow
    sub_sigma{i,1} = subid{idx};
    sub_sigma{(i+1),1} = subid{idx};
    sub_sigma{(i+2),1} = subid{idx};
    idx = idx + 1;
    for j = 1:nsigma
        sub_sigma{(i+j-1),2} = sigmas(j,:);
    end
end

%%% Reassign subid and sigma based on job array counter
% Retrieve SGE_TASK_ID from system (job array index)
batch_idx = getenv('SGE_TASK_ID');

% If this is a job array, then batch_idx will not be empty.
if ~isempty(batch_idx)
    % Convert from ASCII to double
    batch_idx = str2double(batch_idx);
    % Retrieve corresponding row from sub_sigma
    [subid, gsigma] = sub_sigma{batch_idx, :};
% Otherwise, set the Gaussian sigma manually
else
    gsigma = [7, 9, 11];
end


disp({'The subject ID is', subid})
sprintf('The Gaussian sigma array = [%i, %i, %i]', gsigma)

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