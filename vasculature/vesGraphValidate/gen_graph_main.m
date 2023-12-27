%% Main script for creating graph from segmentation
% This script will call the function for generating a graph. This will also
% remove loops from the graph.

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

%% Lookup table of pia pixel intensity reference values
subid = {'AD_10382', 'AD_20832', 'AD_20969',...
             'AD_21354', 'AD_21424',...
             'CTE_6489','CTE_6912',...
             'CTE_7019','CTE_8572','CTE_7126',...
             'NC_6047', 'NC_6839',...
             'NC_6974', 'NC_7597',...
             'NC_8095', 'NC_8653',...
             'NC_21499','NC_301181'};

%% Import files
if ispc
    % Laptop directory structure
    dpath = ['C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\' ...
            'test_data\Ann_Mckee_samples_10T\'];
    % Subfolder containing data
    subdir1 = '\dist_corrected\volume\';
elseif isunix
    % Upper level directory
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T';    
    % Subfolder containing data
    subdir1 = '/dist_corrected/volume/';
end

% Filename with PSOCT scattering tissue volume
vol_name = 'ref_4ds_norm_inv'; 
% Combined segmentations subfolder
subdir2 = 'combined_segs';
% Combined segmentation subfolder
segdir = 'gsigma_3-5-7_5-7-9_7-9-11';
% Filename with combined segmentations
seg_name = 'seg_masked';


%% Iterate through subjects. Generate graph
for ii = 1:length(subid)
    %%% Import combined segmentation file
    fpath = fullfile(dpath, subid{ii}, subdir1, subdir2, segdir, strcat(seg_name,'.mat'));
    tmp = load(fpath,'segm');
    seg = tmp.segm;
    
    %%% Pre-process combined segmentation prior to graphing
    % Remove small segments (fewer than voxmin)
    voxmin = 30;
    seg = rm_short_vessels(seg, voxmin);
    % Remove large islands (likley pia boundary false positives)
    % Some of the segmentations contain large blobs in the pia. I need to
    % write a function to identify + remove these.

    % Smooth segmentation (Gaussian)
    sigma = 3;
    seg = imgaussfilt3(seg,sigma);

    %%% Convert to graph
    % Set the voxel size
    if regexp(vol_name, '4ds')
        vox_dim = [12, 12, 15];
    elseif regexp(vol_name, '10ds')
        vox_dim = [30, 30, 35];
    else
        vox_dim = [30, 30, 35];
    end
    %%% Initialize graph + remove loops + save output
    graph_path = fullfile(dpath, subid{ii}, subdir1, subdir2, segdir);
    % Visualization boolean
    viz = false;
    seg_graph_init(seg, vox_dim, graph_path, seg_name, viz);

end









