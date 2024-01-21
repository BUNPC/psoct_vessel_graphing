%% Apply the mask.tif to the volume.tif
% In some cases, the original volume still contains the background signal.
% In the case that there exists a mask.tif file for this volume, this
% script will call a function to apply the mask to the volume.

clear; clc; close all;

%% Add top-level directory of code repository to path
% Print current working directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
    %%% Init paralell pool
    % Set # threads = # cores for job
    NSLOTS = str2num(getenv('NSLOTS'));
    maxNumCompThreads(NSLOTS);
    % Check to see if we already have a parpool, if not create one with
    % our desired parameters
    poolobj = gcp('nocreate');
    if isempty(poolobj)
	    myCluster=parcluster('local');
	    % Ensure multiple parpool jobs don't overwrite other's temp files
	    myCluster.JobStorageLocation = getenv('TMPDIR');
	    poolobj = parpool(myCluster, NSLOTS);
    end
end
% Remove the two sub folders to reach parent
% (psoct_human_brain\vasculature\vesSegment)
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));

%% Initialize subject ID lists
%%% All subjects to analyze
subid = {'AD_10382', 'AD_20832', 'AD_20969',...
         'AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912',...
         'CTE_7019', 'CTE_8572','CTE_7126',...
         'NC_6839',  'NC_6974', 'NC_8653',...
         'NC_21499', 'NC_301181'};

%% Initialize directories and filenames

%%% Directories 
% Upper level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T';    
% Subfolder with normalized volume
vdir = '/dist_corrected/volume/';
% Subfolder with graph data
gdir = ['/dist_corrected/volume/combined_segs/' ...
    'gsigma_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];

%%% Filenames (.MAT)
gdata = 'seg_refined_masked_graph_data.mat';
gout = 'seg_refined_masked_loopsrm_graph_data.mat';
vname = 'ref_4ds_norm_inv_refined_masked.tif';

%% Iterate subjects, threhsold prob, combine segs, apply mask.
parfor (ii = 1:length(subid), NSLOTS)
    %%% Debugging information
    fprintf('---------Starting Subject %s---------\n',subid{ii})
    
    % Import volume
    fin = fullfile(dpath, subid{ii}, vdir, vname);
    vol = TIFF2MAT(fin);
    % Import graph
    fin = fullfile(dpath, subid{ii}, gdir, gdata);
    g = load(fin);
    d = g.Data;
    g = g.Data.Graph;
    nodes = d.Graph.nodes;
    edges = d.Graph.edges;
    
    %%% Remove loops
    fprintf('removing loops\n')
    % Move to mean minimum voxel intensity
    v_min = 0.99;
    % Down sample search radius
    delta = 6;
    % # iterations for mv2mean function in for-loop iteration in rm_loops
    mv_iter = 1;
    % Variable for displaying debugging windows
    viz = false;
    [node_rm, edges_rm] =...
        rm_loops(nodes, edges, vol, delta, v_min, mv_iter, viz);
    
    %%% Reinitialize graph data after loop removal
    fprintf('reinitializing graph\n')
    g.nodes = node_rm;
    g.edges = edges_rm;
    Data = init_graph(g);

    %%% Save graph with loops removed
    fprintf('saving graph\n')
    fout = fullfile(dpath, subid{ii}, gdir, gout);
    save_graph(fout, Data);    
    
    %%% Debugging Info
    fprintf('---------Finished Subject %s---------\n\n',subid{ii})
end

function save_graph(fout, Data)
% Save volume or segmentation
% INPUTS:
%   fout (string): filepath output
%   vol (matrix): segmentation or volume to save
save(fout, 'Data', '-v7.3');

end
