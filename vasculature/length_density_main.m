%% Vessel length density Analysis
%{
This script is for analyzing the length density in the dataset:
/projectnb/npbssmic/ns/Ann_Mckee_samples_10T

There is a .MAT file named:
ref_4ds_norm_inv_segment_sigma1_thresh0.24_mask40_graph_seglen_um.mat

This .MAT is stored in each of the AD_* and NC_* subjects in the subfolder:
/projectnb/npbssmic/ns/Ann_Mckee_samples_10T/[subject_ID]/dist_corrected/volume/

This .MAT contains an array with the length of each segment in the volume.
These segment lengths are in microns. In order to find the total length of
an entire sample, we just need to sum all elements of this array.

Next, we need to calculate the total volume of each sample. This can be
achieved by finding the non-zero voxels in each slice of the file
"ref_4ds_norm.btf" and then multiplying by the voxel size [12, 12, 15] um.

There are more notes in the code below.

%}


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

%% Data Directory on SCC
% Path to top-level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_10T/';
% Complete subject ID list for Ann_Mckee_samples_10T
subid = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_8572', 'CTE_7126',...
         'NC_21499', 'NC_6047', 'NC_6839', 'NC_6974', 'NC_7597',...
         'NC_8095', 'NC_8653'};
subdir = '/dist_corrected/volume/';
% Volume filename (this will be the same for each subject)
volname = 'ref_4ds.btf';
% Segment length filename 
segname = 'ref_4ds_norm_inv_segment_sigma1_thresh0.24_mask40_graph_seglen_um.mat';

% Voxel dimensions
vox_dim = [12, 12, 15];

%% Loop through subjects
for ii = 1:length(subid)
    %%% Calculate volume of sample from the 
    % Define entire filepath 
    fullpath = fullfile(dpath, subid{ii}, subdir);
    filename = strcat(fullpath, volname);
    % Convert .tif to .MAT
    vol = TIFF2MAT(filename);
    % TODO:
    %   find the non-zero area of each slice. 
    %   multiply this area by the voxel dimensions.
    %   repeat for each slice.
    %   sum volume of each slice to calculate total volume.

    %%% Calculate total length of segments for each sample
    % Define entire filepath 
    fullpath = fullfile(dpath, subid{ii}, subdir);
    filename = strcat(fullpath, segname);
    % Convert .tif to .MAT
    seglen_um = load(filename);
    % TODO:
    %   calculate sum of all elements in seglen_um

    %%% Calculate length density
    % TODO:
    %   Density = (total length) / (total volume)
    %   save this output to a struct for the subject ID

end

%% Generate Histogram of Segment Length Density
% Compare NC vs. AD


