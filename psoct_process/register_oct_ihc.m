%% Register the OCT images with the pathological staining
% Outline:
% - co-register staining (AB & p-taue) with OCT
% - create heat map of staining intensity
% - create heat map of OCT vascular metrics
% - perform regression analysis between staining & vascular metric

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
    ptau_path = '/projectnb/npbssmic/pantinew/tau_amyloid_images/AT8/';
    ab_path = '/projectnb/npbssmic/pantinew/tau_amyloid_images/Ab/';
    % Subfolder of OCT images
    refdir = '/dist_corrected/volume/';
    % Subfolder of graph data
    subdir = ['/dist_corrected/volume/combined_segs/' ...
        'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];
end



% Filename to parse (this will be the same for each subject)
seg_name = 'seg_refined_masked.tif';
graph_name = 'seg_refined_masked_graph_data.mat';

%%% Complete subject ID list for Ann_Mckee_samples_10T
% subid = {'AD_10382', 'AD_20832', 'AD_20969',...
%          'AD_21354', 'AD_21424',...
%          'CTE_6489','CTE_6912',...
%          'CTE_7019','CTE_8572','CTE_7126',...
%          'NC_6047', 'NC_6839',...
%          'NC_6974', 'NC_7597',...
%          'NC_8095', 'NC_8653',...
%          'NC_21499','NC_301181'};
subid = {'CTE_6489'};

%%% Staining slice number
pathology = struct();
pathology.(subid{1}).ab = 14;
pathology.(subid{1}).ptau = 8;

%% Load/calculate graph metrics

% Define entire filepath to graph
fullpath = fullfile(dpath, subid, subdir);

% Load graph
gmet = load(char(fullfile(fullpath, graph_name)));
gmet = gmet.Data.Graph;
nodes = gmet.nodes;
edges = gmet.edges;

%%% Branch Density
% Calculate degree of nodes: 1 = end. 2 = middle. >2 = branch
g = graph(edges(:,1), edges(:,2));
d = degree(g);
% Identify branch points (degree > 2)
branch_node = nodes(d>2,:);
% TODO: calculate regional branch density w/ coordinates of branch nodes

%%% Fraction Volume


%%% Length Density


%%% Tortuosity



%% Load the IHC staining (AT8 for p-tau)
% Different slices were stained for either AB or p-tau. In each case, the
% brightfield image (with all color channels) shows the edges of the
% tissue. These should be used to register the staining to the OCT images.

%%% Load brightfield image of tissue
bfield = char(fullfile(ptau_path,'/CTE 6489 AT8/CTE_6489_brightfield_20ds.tif'));
bfield = Tiff(bfield, 'r');
bfield = read(bfield);
figure; imshow(bfield);

%%% Load p-tau
ptau = char(fullfile(ptau_path,'/CTE 6489 AT8/CTE_6489_AT8_20ds.tif'));
ptau = Tiff(ptau, 'r');
ptau = read(ptau);
% Convert p-tau to grayscale
ptau_bw = rgb2gray(ptau);
figure; imshow(ptau);
figure; imagesc(ptau_bw); colorbar; title('p-tau w/ background')
% Remove background (purely white pixels (255)
ptau_bw(ptau_bw==255) = 0;
figure; imagesc(ptau_bw); colorbar; title('p-tau w/o background')

%%% Load the respective OCT slice
slicen = pathology.(subid{1}).ptau;
refname = strcat('Ref_BASIC',num2str(slicen),'.tif');
oct = fullfile(dpath, subid{1}, refdir, refname);
oct = Tiff(oct, 'r');
oct = read(oct);

% Register brightfield and OCT


% Apply registration coordinates to p-tau and OCT


%% Load the IHC staining (Amyloid Beta)
%%% AB Stain
% Load brightfield image of tissue

% Load AB stain channel

% Load the respective OCT slice

% Register brightfield and OCT

% Apply registration coordinates to AB and OCT


%% Heat map - IHC

% Downsample?

% Display 


%% Heat map - OCT metrics (fraction volume)



%% Correlate heat maps - moving window regression





















