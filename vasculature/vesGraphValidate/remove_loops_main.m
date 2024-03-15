%% Main script for removing loops from preexisting graphs
% This script will call rm_loops_parallel for removing loops from graphs.
clear; clc; close all;
dbstop if error
%% Add top-level directory of code repository to path
% Print current working directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Remove the two sub folders to reach parent
% (psoct_human_brain\vasculature\vesSegment)
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));

%% Initialize subject list and directory structure

% Subset of subjects with good segmentation
subid ={'AD_10382', 'AD_20832', 'AD_20969','AD_21354','AD_21424',...
       'CTE_6489','CTE_6912','CTE_7019','CTE_7126',...
       'NC_6839','NC_6974','NC_8653','NC_21499','NC_301181'};
% Top-level Directory Path
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
% Subfolder with OCT volume
subdir = '/dist_corrected/volume/';
vdata = 'ref_4ds_norm_inv_refined_masked.tif';
% Sub-subfolder with OCT volume segmentation and graph
sigdir = '/combined_segs/gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p20/';
seg_name = 'seg_refined_masked.tif';
gdata = 'seg_refined_masked_graph_data.mat';

%% Iterate through subjects. Generate graph
%%% Import combined segmentation file
for ii = 8:length(subid)
    %%% Load Graph
    % Subject ID in loop
    sub = subid{ii};
    sprintf('\n\nStarting subject %s\n\n',sub)
    % Load graph data
    Data = load(fullfile(dpath, sub, subdir, sigdir, gdata), 'Data');
    Data = Data.Data;
    nodes = Data.Graph.nodes;
    edges = Data.Graph.edges;
    
    %%% Load volumetric information
    vol = TIFF2MAT(fullfile(dpath, sub, subdir, vdata));
    
    %%% Variables for Loop Removal
    % Down sample search radius
    delta = 4;
    % Move to mean minimum voxel intensity
    v_min = 0.99;
    % # iterations for mv2mean function in for-loop iteration in rm_loops
    mv_iter = 1;
    % Boolean for visualizing debugging graphs
    viz = false;
    
    %%% Call function to remove loops
    [nodes_rm, edges_rm] =...
        rm_loops_parallel(nodes, edges, vol, delta,...
        v_min, mv_iter, viz);
    % Update the graph
    Data.Graph.nodes = nodes_rm;
    Data.Graph.edges = edges_rm;
    % Reinitialize the graph metadata after loop removal
    Data = init_graph(Data.Graph);

    %%% Save output to same folder
    gdata_out = strcat(gdata(1:end-4),'_loops_rm.mat');
    graph_filepath = fullfile(dpath, sub, subdir, sigdir, gdata_out);
    save(graph_filepath,'Data');
    sprintf('\n\nFinished subject %s\n\n',sub)
end