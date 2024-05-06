%% Divide volume into grid and calculate metrics for each cube
% The purpose of this is to subdivide the entire tissue volume in cubes of
% dimensions [M, M, M] cubic millimeters and then generate the graph for
% each cube. The following step is to generate the heat map from these
% cubes.

%% Add top-level directory of code repository to path
clear; clc; close all;
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

%% Initialize paralell pool
% Set # threads = # cores for job
NSLOTS = str2num(getenv('NSLOTS'));
maxNumCompThreads(NSLOTS);
% Check to see if we already have a parpool, if not create one with
% our desired parameters
% poolobj = gcp('nocreate');
% if isempty(poolobj)
%     myCluster=parcluster('local');
%     % Ensure multiple parpool jobs don't overwrite other's temp files
%     myCluster.JobStorageLocation = getenv('TMPDIR');
%     poolobj = parpool(myCluster, NSLOTS);
% end

%% Initialize subject ID lists
%%% All subjects to analyze
subid = {'AD_10382', 'AD_20832', 'AD_20969','AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912','CTE_7019','CTE_7126',...
         'NC_6839',  'NC_6974', 'NC_8653','NC_21499', 'NC_301181'};

%%% stacks that require truncating
% Last depth to retain for each stack
zmins = [187, 165, 242, 165, 220,...
        242, 110, 198, 198,...
        176, 198, 165, 165, 220];
% Create dictionary to store last image in stack
d = dictionary(subid, zmins);

%% Initialize directories, filenames, parameters

%%% Directories 
% Upper level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T';    
% Subfolder with normalized volume
subdir = '/dist_corrected/volume';
% Subfolder containing non-normalized ref files
subdir1 = '/dist_corrected/volume/ref';
% Combined segmentation subfolder
segdir = '/combined_segs/gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/';
% Mask subfolder
mdir = '/dist_corrected/volume/ref/masks';
% Metrics output path
mpath = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
    'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];

%%% Filenames
% Masked combined segmentations (.MAT)
seg_name = 'seg_refined_masked.mat';

%%% Subvolume parameters
% Size of each cube (cubic microns)
cube_vol_vox = 1e9;
% Size of each voxel (microns)
vox = [12, 12, 15];
% Compute number of voxels in x,y,z for each cube
cube_side = nthroot(cube_vol_vox,3);
n_x = floor(cube_side ./ vox(1));
n_y = floor(cube_side ./ vox(2));
n_z = floor(cube_side ./ vox(3));

%% Generate grid for each subject
% TODO:
% - Calculate branch density (from entire graph)
% - Overlay tissue mask on imagesc of metrics

%%% Struct for storing heat map
heatmap = struct();

%%% Iterate over each subject
for ii = 1:length(subid)
    %% Generate Subvolumes
    % Load masked segmentation
    masked_seg = fullfile(dpath,subid{ii},subdir,segdir,seg_name);
    seg = load(masked_seg);
    seg = seg.seg;  

    % Calculate number of cubes in x, y, z
    Nx = ceil(size(seg,1) ./ n_x);
    Ny = ceil(size(seg,2) ./ n_y);
    Nz = ceil(size(seg,3) ./ n_z);
    Ntot = Nx * Ny * Nz;

    % Find remainder for each dimension
    rem_x = rem(size(seg,1), n_x);
    rem_y = rem(size(seg,2), n_y);
    rem_z = rem(size(seg,3), n_z);

    %% Iterate over the segmentation
    % Initialize the volume fraction heat map matrix 
    vf_mat = zeros(size(seg,1), size(seg,2), Nz);
    ld_mat = zeros(size(seg,1), size(seg,2), Nz);
    % Calculate the size of each cube (in voxels)
    cube_vol_vox = n_x * n_y * n_z;
    cube_vol_um = cube_vol_vox * vox(1) * vox(2) * vox(3);
    
    % Iterate over the z-axis
    for z = 1:n_z:size(seg,3)
        % Iterate over rows
        for x = 1:n_x:size(seg,1)
            % Iterate over columns
            for y = 1:n_y:size(seg,2)
                %% Crop segmentation into cube
                % Initialize end indices for each axis
                xf = x + n_x - 1;
                yf = y + n_y - 1;
                zf = z + n_z - 1;
                % Take minimum of matrix dimensions and end indices
                xf = min(xf, size(seg,1));
                yf = min(yf, size(seg,2));
                zf = min(zf, size(seg,3));
                % Take cube from segmentation
                cube = seg((x:xf), (y:yf), (z:zf));
                
                %% Volume Fraction
                % Calculate volume fraction of cube
                vf = sum(cube(:)) ./ cube_vol_vox;
                % Assign to the vf matrix
                vf_mat((x:xf), (y:yf), (z:zf)) = vf;
                
                %% Calculate length density if cube has segmentation
                % Check for segmentation within cube
                if sum(cube(:)) > 1
                    %%% Create graph with subvolume (cube)
                    [graph, ~] = seg_to_graph(cube, vox);
                    % Move to mean minimum voxel intensity
                    v_min = 0.99;
                    % initial search radius in down sample function (voxels)
                    delta = 6;
                    % # iterations for mv2mean function
                    mv_iter = 1;
                    % Visualize debugging
                    viz = false;
                    % Call remove loops
                    [nodes_rm, edges_rm] =...
                        rm_loops_parallel(graph.nodes, graph.edges, cube,...
                                          delta, v_min, mv_iter, viz);
                    % Update graph with nodes and edges
                    graph.nodes = nodes_rm;
                    graph.edges = edges_rm;
                    
                    %%% If there is only 1 edge, then skip metadata. There
                    % is a bug in the graphing code that raises an error
                    % when there is only one edge. This is legacy code and
                    % difficult to debug. Instead, just calculate Euclidean
                    % distance of the lone edge.
                    if size(edges_rm,1) > 1
                        % Compute the graph metadata
                        [Data] = init_graph(graph);
        
                        %%% Calculate length density with subvolume (cube)
                        seglen_um = Data.Graph.segInfo.segLen_um;
                        seglen_tot_um = sum(seglen_um);
                        ld = seglen_tot_um ./ cube_vol_um;
                    elseif size(edges_rm,1) == 1
                        %%% Euclidean Distance
                        n1 = nodes_rm(edges_rm(1),:) .* vox;
                        n2 = nodes_rm(edges_rm(2),:) .* vox;
                        d = sqrt( (n1(1) - n2(1)).^2 +...
                                  (n1(2) - n2(2)).^2 +...
                                  (n1(3) - n2(3)).^2);
                        ld = d ./ cube_vol_um;
                    else
                        %%% in this case there is no graph
                        ld = 0;
                    end
                    %%% Add length density to matrix
                    ld_mat((x:xf), (y:yf), (z:zf)) = ld;
                end
            end
        end
    end
   
    %% Separate LD and VF matrices into z-dimension slices by cube size
    % Create array of z-dimension indices to retain
    idx_z = 1 : n_z : size(seg,3);
    
    % Initialize matrices to store volume fraction and length density
    heatmap_vf = zeros(size(vf_mat,1), size(vf_mat,2), length(idx_z));
    heatmap_ld = zeros(size(ld_mat,1), size(ld_mat,2), length(idx_z));
        
    % Keep just the idx_z frames
    for j = 1:length(idx_z)
        % Extract volume fraction
        heatmap_vf(:,:,j) = vf_mat(:,:,idx_z(j));
        % Extract length density
        heatmap_ld(:,:,j) = ld_mat(:,:,idx_z(j));
    end

    % Add volume fraction and length density to struct
    sub = subid{ii};
    heatmap.(sub).vf = heatmap_vf;
    heatmap.(sub).ld = heatmap_ld;


    %% Older method for segmenting into cube
    %{
    %%% Calculate number of cubes in x, y, z
    rem_x = rem(size(seg,1), n_x);
    rem_y = rem(size(seg,2), n_y);
    rem_z = rem(size(seg,3), n_z);

    %%% Pad segmentation in X,Y to enable whole number of cubes
    pad_x = n_x - rem_x;
    pad_y = n_y - rem_y;
    seg = cat(1, seg, zeros(pad_x, size(seg,2), size(seg,3)));
    seg = cat(2, seg, zeros(size(seg,1), pad_y, size(seg,3)));
    
    %%% Truncate the z-dimension to enable whole number of cubes
    % Find number of cubes in z dimension
    N_z = floor(size(seg,3) ./ n_z);
    z_idx = N_z * n_z;
    seg_cubes = seg(:,:,1:z_idx);
    % Place remaining z-dimension into separate matrix
    seg_cubes_rem = seg(:,:,(z_idx+1):end);

    %%% Divide into cubes of [n_x, n_y, n_z] voxels
    seg_cubes = reshape(logical(seg_cubes), n_x, n_y, n_z, []);
    seg_cubes_rem = reshape(logical(seg_cubes_rem), n_x, n_y, rem_z, []);

    %% Calculate volume fraction for each subvolume
    
    % Initialize matrix to store volume fraction
    x_subvol = size(seg_cubes,1);
    y_subvol = size(seg_cubes,2);
    z_subvol = size(seg_cubes,3);
    n_subvol = size(seg_cubes,4);
    vf_subvol = zeros(x_subvol, y_subvol, z_subvol, n_subvol);

    % Calculate the size of each cube (in voxels)
    cube_vol = size(seg_cubes,1) .* size(seg_cubes,2) .* size(seg_cubes,3);

    % Iterate over each cube
    for j = 1:size(seg_cubes,4)
        % Take j_th cube
        cube = seg_cubes(:,:,:,j);
        % Calculate volume fraction of cube
        vf = sum(cube(:)) ./ cube_vol;
        % Add volume fraction to the matrix
        vf_subvol(:,:,:,j) = vf;
    end

    % Reshape cubes back into matrix
    dim = [size(seg,1), size(seg,2), z_idx];
    vf_subvol = reshape(vf_subvol, dim);
    % Truncate to remove zero padding
    vf_subvol = vf_subvol(1:(end-pad_x), 1:(end-pad_y), :);
    %}

end

% Save the heatmap
heat_out = fullfile(mpath, 'heatmap.mat');
save(heat_out,'heatmap','-v7.3');

%% Plot and save the heat maps
function plot_save_heatmap(Nframes, heatmap, tstr, dpath, fname)
% PLOT_SAVE_HEATMAP use imagesc and set background = 0
% INPUT
%   Nframes (int): number of frames in z dimension
%   heatmap (double): heatmap of vascular metric
%   tstr (string): figure title
%   dpath (string): data directory path
%   fname (string): name of figure to save

%%% Iterate over frames in z dimension
for j = 1:length(Nframes)
    % Take the j_th frame from the length density
    ld = heatmap(:,:,j);
    cmap_min = min(ld(ld(:)>0));
    cmap_max = max(ld(:));
    imagesc(ld);
    cmap = jet(256);
    colormap(cmap);
    clim(gca, [cmap_min, cmap_max]);
    % Make value 0 black
    cmap(1,:) = [0,0,0];
    colormap(cmap);
    colorbar;
    
    %%% Save figure as PNG
    fout = fullfile(dpath, fname);
    saveas(gca, fout,'png');

end
end