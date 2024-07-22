%% Remove loops main script
% Purpose: remove the loops from the graphs of the entire tissue volume for
% each CAA subject

clear; close all; clc;

%% Add top-level directory of code repository to path
% Start in current directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Truncate path to reach top-level directory (psoct_vessel_graphing)
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));

%% Initialize data path for CAA graphs
% Path to top-level directory
dpath = '/projectnb/npbssmic/ns/CAA/';
% Subfolder containing data
subdir = 'segmentations/';
% Filename to parse (this will be the same for each subject)
fnames = {'caa22-occipital_vessels-masked',...
          'caa25-occipital_vessels-masked',...
          'caa26-frontal_vessels-masked',...
          'caa26-occipital_vessels-masked'};
% subjects with loops remaining in graph
subids = {'caa22/occipital/','caa25/occipital/',...
         'caa26/frontal/', 'caa26/occipital/'};

%%% Select the subject ID and filename based on the batch index
% Retrieve SGE_TASK_ID from system (job array index)
batch_idx = getenv('SGE_TASK_ID');
% Debugging: manually set the index
if strcmp(batch_idx,'undefined')
    batch_idx = '4';
end
% Convert from ASCII to double
batch_idx = str2double(batch_idx);
% Retrieve corresponding row from sub_sigma
subid = subids{batch_idx};
fname = fnames{batch_idx};

%% Import the graph and segmentation

% Path to graph
dpath = fullfile(dpath, subid, subdir);
Data = append(fname,'_graph_data.mat');
Data = strcat(dpath, Data);

% Load graph
Data = load(Data);
Data = Data.Data;
graph = Data.Graph;
nodes = graph.nodes;
edges = graph.edges;
vox = graph.vox;

% Load segmentation
seg = Data.angio;

%% Remove loops from graph

%%% Input parameters for Remove Loops Parallel function 
% Move to mean minimum voxel intensity
v_min = 0.99;
% Down sample search radius
delta = 4;
% # iterations for mv2mean function in for-loop iteration in rm_loops
mv_iter = 1;
% Boolean for visualizing debugging graphs
viz = false;

%%% Call the function to remove loops
[nodes_rm, edges_rm] =...
    rm_loops_parallel(nodes, edges, seg, delta, v_min, mv_iter, viz);
% Reassign the new nodes and edges
graph.nodes = nodes_rm;
graph.edges = edges_rm;

%% Save the output to the same directory as input data
% Define the output filename
fout = append(fname,'_rmloop_graph_data.mat');
fout = fullfile(dpath, fout);
% Save output as struct
save(fout, 'graph', '-v7.3');