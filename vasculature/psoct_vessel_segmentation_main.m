%% Main file for calling segmentation functions
% Author: Mack Hyman
% Date Created: March 16, 2023
%
% Detailed Description
%{
This script performs the following:
- segment the original volume
- apply a mask to the segmentation
- convert segmentation to graph
To Do:
- find optimal range for remove_mask_islands
- prune graph (remove loops and unterminated segments)
    - remove loops ()
    - remove segments ()
%}
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
    dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\';
    % Subject IDs
    subid = {'CTE_7019'};
    subdir = '\dist_corrected\volume\';
    % Filename to parse (this is test data)
    fname = 'ref_4ds_norm_inv_cropped';
    % filename extension
    ext = '.tif';
    % sigma for Gaussian smoothing
    gsigma = [7, 9, 11];

%%% Computing cluster (SCC)
elseif isunix
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

    
end

%% Initialization parameters (same for both 

%%% Assign PS-OCT voxel dimension [x, y, z] according to downsample factor
% Downasample factor = 4 --> Voxel = [12, 12, 15] micron
% Downasample factor = 10 --> Voxel = [30, 30, 35] micron
% 2P microscopy pixel will always be [2, 2] micron
if regexp(fname, '4ds')
    vox_dim = [12, 12, 15];
elseif regexp(fname, '10ds')
    vox_dim = [30, 30, 35];
else
    vox_dim = [30, 30, 35];
end

%%% Size of the Gaussian kernel. This should be a 3-element array of
% positive, odd integers. Default size is 2*ceil(2*gsigma)+1
gsize = 2.*ceil(2.*gsigma)+1;

%%% Minimum fringi filter probability to classify voxel as vessel
min_prob = 0.20:0.01:0.26;

%%% A segment with < "min_conn" voxels will be removed
min_conn = 30;

%%% Array (or single value) of radii for eroding the mask
radii = 40;

%%% Boolean for converting segment to graph (0 = do not convert. 1 = convert)
graph_boolean = 1;

for ii = 1:length(subid)
    %% Load raw volume (TIF) and convert to MAT
    % Define entire filepath 
    fullpath = fullfile(dpath, subid{ii}, subdir);
    filename = strcat(fullpath, strcat(fname, ext));
    % Convert .tif to .MAT
    vol = TIFF2MAT(filename);

    %%% Create subfolder for Gaussian sigma and kernel size
    % Create string of Gaussian sigmas
    gsigma_str = num2str(gsigma);
    % Replace spaces with hyphens
    gsigma_str = strrep(gsigma_str, '  ', '-');
    gsigma_subfolder = strcat('gsigma_',gsigma_str);
    
    % Create string of Gaussian kernel sizes
    gsize_str = num2str(gsize);
    % Replace spaces with hyphens
    gsize_str = strrep(gsize_str, '  ', '-');
    gsize_subfolder = strcat('_gsize_',gsize_str);

    % concatenate sigma and kernel into single directory
    subfolder = strcat(gsigma_subfolder, gsize_subfolder);

    % Create string for entire directory path to subfolder
    fullpath = fullfile(fullpath, subfolder);
    
    % Create subfolder with Gaussian sigma and kernel size
    if ~exist(fullpath, 'dir')
       mkdir(fullpath)
       % Add metadata text file
       metadata = {strcat('gaussian sigma =  ', num2str(gsigma)),...
           strcat('gaussian kernel =  ', num2str(gsize)),...
           strcat('minimum probability =  ', num2str(min_prob)),...
           };
       writelines(metadata, fullfile(fullpath, 'metadata.txt'));
    end


    %% Segment volume
    % convert volume to double matrix
    vol = double(vol);
    % Segment volume. Threshold with first element of probability matrix.
    [pmat, seg] = vesSegment(vol, gsigma, gsize, min_prob(1), min_conn);
    % Save probability map for posterity
    fout = strcat(fullfile(fullpath, 'probability_map'), '.mat');
    save(fout, 'pmat', '-v7.3');
    
    for j = 1:length(min_prob)
        %%% Threshold probability matrix with min_prob array
        I_seg = pmat;
        I_seg(pmat < min_prob(j)) = 0;
        I_seg(pmat >= min_prob(j)) = 1;
        % Convert binary matrix to unsigned 8-bit to save memory
        I_seg = uint8(I_seg);
        % Remove segments with fewer than voxmin connected voxels
        I_seg = rm_short_vessels(I_seg, min_conn);

        %%% Save unmasked & thresholded segmentation to TIF
        % Create filename for probability
        fname_seg = strcat(fname,'_segment_pmin_',num2str(min_prob(j)));
        fout = strcat(fullfile(fullpath, fname_seg), '.tif');
        segmat2tif(I_seg, fout);

        %%% Overlay mask volume (grayscale) and unmasked segmentation (green)
        % Create output filename
        overlay_name = strcat(fname_seg, '_overlay.tif');
        overlay_fout = fullfile(fullpath, overlay_name);
        % Call function to overlay mask and segmentation
        overlay_vol_seg(vol, I_seg, 'green', overlay_fout);

        if graph_boolean
            seg_graph_init(I_seg, vox_dim, fullpath, fname_seg);
        end

        %% Mask segmented volume (remove erroneous vessels) & Convert to Graph
        % The function for creating the mask requires a radius. This for-loop will
        % iterate over an array of radii. For each radius, it will create a mask,
        % apply the mask to the segmentation volume, and save the output.
        % If the graph_boolean is true (1), then the masked segmentation will be
        % converted to a graph.
    
        % Create 3D mask from original volume
        mask = logical(vol);
        for k = 1:length(radii)
            %%% Apply mask and save .MAT and .TIF
            [I_seg_masked] = mask_segments(I_seg, mask, radii(k),...
                                            fullpath, fname_seg);
            
            %%% Overlay mask volume (grayscale) and segmentation (green)
            % Create output filename
            overlay_name = strcat(fname_seg,'_mask',num2str(radii(k)),'_overlay.tif');
            overlay_fout = fullfile(fullpath, overlay_name);
            % Call function to overlay mask and segmentation
            overlay_vol_seg(vol, I_seg_masked, 'green', overlay_fout);

            %%% Convert masked segmentation to graph
            if graph_boolean
                fname_masked = strcat(fname_seg, '_mask_', num2str(radii(k)));
                seg_graph_init(I_seg_masked, vox_dim, fullpath, fname_masked);
            end
        end
    end
end

function seg_graph_init(seg, vox_dim, fullpath, fname_seg)
% Initialize graph from segmentation
% INPUTS:
%   seg (mat): segmentation matrix
%   vox_dim (array): 3-element array of voxel dimensions (microns)
%

%%% Convert segmentation to graph (just the nodes and segments)
graph_nodes_segs = seg_to_graph(seg, vox_dim);

%%% Initialize graph metadata (Graph.Data)
[Data] = init_graph(graph_nodes_segs);

%%% Append "angio" data (segmentation matrix)
% Rearrange [x,y,z] --> [z,x,y]. This is the standard in
% the graph validation GUI.
angio = permute(seg, [3,1,2]);
Data.angio = angio;

% Create new filename for graph and add .MAT extension
fname_graph = strcat(fname_seg, '_graph_data.mat');
fout = fullfile(fullpath, fname_graph);
save(fout,'Data', '-v7.3');
end

%% Apply Mask
function [I_seg_masked] = mask_segments(I_seg, mask, radius, fullpath, fname)
% Remove the edges labeled as vessels.
%   INPUTS:
%       I_seg (matrix) - output of segmentation function
%       mask (matrix) - unsegmented volume converted to logicals
%       radius (double array) - radius of disk for eroding the mask
%       fullpath (string) - absolute directory for saving processed data
%       fname (string) - filename prior to applying mask
%   OUTPUTS:
%       I_seg_masked (matrix) - I_seg with boundaries eroded to remove
%           erroneously labeled vessels.

%%% Erode mask to remove small pixels on border that are not part of volume
se = strel('disk', radius);
mask = imerode(mask, se);

%%% Remove islands of pixels from mask
% Range of object size to keep
range = [1e4, 1e8];
mask = remove_mask_islands(mask, range);

%%% Apply mask to segmentation volume
% Convert from logical back to uint8 for matrix multiplication
mask = uint8(mask);
% Element-wise multiply mask and volume
I_seg_masked = apply_mask(I_seg, mask);
% Convert segmentation output to uint8
I_seg_masked = uint8(I_seg_masked);

%%% Save segmented/masked volume as .MAT and .TIF
% Convert masked image back to tif
tmp_fname = strcat(fname,'_mask', num2str(radius));
fout = fullfile(fullpath, tmp_fname);
fout = strcat(fout, '.tif');
segmat2tif(I_seg_masked, fout);
% Save vessel segment stack as .MAT for the next step (graph recon)
fout = fullfile(fullpath, tmp_fname);
fout = strcat(fout, '.mat');
save(fout, 'I_seg_masked', '-v7.3');

end