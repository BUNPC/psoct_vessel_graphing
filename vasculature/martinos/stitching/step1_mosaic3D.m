%% Step 1 of pre-processing
%{
This script will stitch each slice.
For each data collection run:
    - change the ParameterFile for each "process_run#"
    - execute this script
To Do:
    - parameterize function call for each run with for-loop
%}

%% Initial Parameters
clear; clc;
%%% Add top-level directory to path
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Truncate path to reach top-level directory (psoct_vessel_graphing)
topdir = mydir(1:idcs(end-2));
addpath(genpath(topdir));

% Specify imaging modality
modality = 'mus';
% Path to parameter file
ParameterFile = '/autofs/cluster/octdata2/users/Chao/caa/caa_6/frontal/process_run1/Parameters.mat';
stitch_xy(ParameterFile, modality);

%%% WIP code for serializing function call
% ParameterFile = {'', '', ''};
% nruns = length(ParameterFile);
% for run = 1:nruns
%     stitch_xy(ParameterFile{run}, modality);
% end

