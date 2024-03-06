%% Combine the b-scans (ref#.mat) to remake the c-scan
% Each subdirectory ".../[subjectID]/dist_corrected/volume/ contains the
% b-scans in the format ref#.mat. It is necessary to remake the c-scans to
% make a mask for removing the tissue boundary. This is necessary for
% removing false positives along the border of the tissue volumes.

clear; clc; close all;

%% Add top-level directory of code repository to path
% This allows Matlab to find the functions in the project folders

% Start in current directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Truncate path to reach top-level directory (psoct_vessel_graphing)
topdir = mydir(1:idcs(end));
addpath(genpath(topdir));

%% Initialize data path for linux or personal machine (debugging)

%%% Local machine
if ispc
    dpath = ['C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures'...
            '\test_data\Ann_Mckee_samples_10T\'];
%%% Computing cluster (SCC)
elseif isunix
    % Path to top-level directory
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
end

% Subfolder containing data
subdir = '/dist_corrected/volume/';
% Filename to parse (this will be the same for each subject)
fname = 'ref';
%%% Complete subject ID list for Ann_Mckee_samples_10T
subid = {'AD_10382', 'AD_20832', 'AD_20969',...
         'AD_21354', 'AD_21424',...
         'CTE_6489','CTE_6912',...
         'CTE_7019','CTE_8572','CTE_7126',...
         'NC_6047', 'NC_6839',...
         'NC_6974', 'NC_7597',...
         'NC_8095', 'NC_8653',...
         'NC_21499','NC_301181'};

%% Load raw volume (TIF) and convert to MAT
% Define entire filepath 
fullpath = fullfile(dpath, subid, subdir);

% Find files with numbers
fpath = fullpath{1};
list = ls(fpath);

% Regular Expression to find mat files with numbers
[~, reindex] = sort( str2double( regexp(list,'\d+', 'match', 'once' )));
% [~, reindex] = sort( str2double( regexp( {A.name}, '\d+', 'match', 'once' )));
A = A(reindex) ;

%
filename = strcat(fullpath, strcat(fname, ext));

%% Find all files with "ref#.mat" in the subfolder