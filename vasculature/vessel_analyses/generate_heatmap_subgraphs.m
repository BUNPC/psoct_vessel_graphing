%% Divide volume into grid and calculate metrics for each cube
% The purpose of this is to subdivide the entire tissue volume in cubes of
% dimensions [M, M, M] cubic millimeters and then generate the graph for
% each cube. The following step is to generate the heat map from these
% cubes.

% TODO:
%{
- modularize the heatmap to operate at a certain depth
- define the respective OCT depth that aligns with the pathology depth
- generate heatmaps and masks for each OCT/pathology depth
%}

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
% Set maximum number of threads equal to number of threads for script
NSLOTS = str2num(getenv('NSLOTS'));
maxNumCompThreads(NSLOTS);

%% Initialize directories, filenames, parameters

%%% All subjects to analyze for vasculature
subid = {'AD_10382', 'AD_20832', 'AD_20969','AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912','CTE_7019','CTE_7126',...
         'NC_6839',  'NC_6974', 'NC_8653','NC_21499', 'NC_8095'};
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
% Graph (without loops) from masked segmentation
graph_name = 'seg_refined_masked_rmloop_graph_data.mat';

%%% Subvolume parameters
% Isotropic cube length (microns)
cube_side = 204;
% Size of each voxel (microns)
vox = [12, 12, 15];
% Compute number of voxels in x,y,z for each cube
n_x = floor(cube_side ./ vox(1));
n_y = floor(cube_side ./ vox(2));
n_z = floor(cube_side ./ vox(3));
% Calculate the size of each cube (in voxels)
cube_vol_vox = n_x * n_y * n_z;
% Calculate the size of each cube (in cubic microns)
cube_vol_um = cube_vol_vox * vox(1) * vox(2) * vox(3);

%%% Struct for storing vascular heat map
heatmap = struct();

%%% Whether to plot non-normalized heatmaps
viz_individual = false;

%% Pathology heatmap: A-beta and p-tau
%{
% Initialize maximum/minimum of all heat maps
ab_max = 0;
ab_min = 1;
pt_max = 0;
pt_min = 1;
% Initialize struct to store pathology heatmap
path_heatmap = struct();
% Title above figures
tstr = {'A-beta','P-tau'};
% Types of stains
stains = {'ab','pt'};
% Filenames of each stain
stain_fname = struct();
stain_fname.AD_10382.ab = 'AD_10382_slice_8_Ab';
stain_fname.AD_10382.pt = 'AD_10382_slice_14_AT8';
stain_fname.AD_20832.ab = 'AD_20832_slice_8_Ab';
stain_fname.AD_20832.pt = 'AD_20832_slice_14_AT8';
stain_fname.AD_21354.ab = 'AD_21354_slice_8_Ab';
stain_fname.AD_21354.pt = 'AD_21354_slice_2_AT8';
stain_fname.AD_21424.ab = 'AD_21424_slice_14_Ab';
stain_fname.AD_21424.pt = 'AD_21424_slice_20_AT8';
stain_fname.CTE_6489.ab = 'CTE_6489_slice_14_Ab';
stain_fname.CTE_6489.pt = 'CTE_6489_slice_8_AT8';
stain_fname.CTE_6912.ab = 'CTE_6912_slice_8_Ab';
stain_fname.CTE_6912.pt = 'CTE_6912_slice_14_AT8';
stain_fname.CTE_7019.ab = 'CTE_7019_slice_8_Ab';
stain_fname.CTE_7019.pt = 'CTE_7019_slice_14_AT8';
stain_fname.CTE_7126.ab = 'CTE_7126_slice_14_Ab';
stain_fname.CTE_7126.pt = 'CTE_7126_slice_8_AT8';
stain_fname.NC_21499.ab = 'NC_21499_slice_14_Ab';
stain_fname.NC_21499.pt = 'NC_21499_slice_8_AT8';
stain_fname.NC_8095.ab = 'NC_8095_slice_8_Ab';
stain_fname.NC_8095.pt = 'NC_8095_slice_2_AT8';
% Minimum thresholds for segmenting plaques (after taking compliment)
stain_thresh = struct();
stain_thresh.AD_10382.ab = 0.2;
stain_thresh.AD_10382.pt = 0.2;
stain_thresh.AD_20832.ab = 0;
stain_thresh.AD_20832.pt = 0.10;
stain_thresh.AD_21354.ab = 0.23;
stain_thresh.AD_21354.pt = 0.11;
stain_thresh.AD_21424.ab = 0.15;
stain_thresh.AD_21424.pt = 0.15;
stain_thresh.CTE_6489.ab = 0.16;
stain_thresh.CTE_6489.pt = 0.13;
stain_thresh.CTE_6912.ab = 0.10;
stain_thresh.CTE_6912.pt = 0.13;
stain_thresh.CTE_7019.ab = 0.15;
stain_thresh.CTE_7019.pt = 0.15;
stain_thresh.CTE_7126.ab = 0.15;
stain_thresh.CTE_7126.pt = 0.15;
stain_thresh.NC_21499.ab = 0.10;
stain_thresh.NC_21499.pt = 0.10;
stain_thresh.NC_8095.ab = 0.10;
stain_thresh.NC_8095.pt = 0.10;
% Subject ID list containing staining
stain_subid = fieldnames(stain_fname);

for ii = 1:length(fields(stain_fname))
    for j = 1:length(stains)
        %%% A-Beta: load pathology, apply mask, calculate resolution
        fpath = fullfile(dpath, stain_subid{ii}, 'stain');
        % Filename of staining
        stain_name = stain_fname.(stain_subid{ii}).(stains{j});
        % Generate mask file name
        mask_name = append(stain_fname.(stain_subid{ii}).(stains{j}),'_mask');
        % Import the stain, mask, and calculate pxiel resolution
        [stain, mask, res] =...
            import_pathology(fpath, stain_name, mask_name);
        
        %%% generate heatmap from A-Beta staining
        th = stain_thresh.(stain_subid{ii}).(stains{j});
        [hm] = pathology_heatmap(res, cube_side, stain, th, mask);
        % Calculate maximum value of all heatmaps
        if strcmp(stains{j},'ab')
            ab_max = max(ab_max, max(hm(mask)));
            ab_min = min(ab_min, min(hm(mask)));
        else
            pt_max = max(pt_max, max(hm(mask)));
            pt_min = min(pt_min, min(hm(mask)));
        end

        %%% Add heatmap of A-beta and p-tau to struct
        path_heatmap.(stain_subid{ii}).(stains{j}).heatmap = hm;
        path_heatmap.(stain_subid{ii}).(stains{j}).mask = mask;
    end
end

%%% Plot and save the pathology heatmaps
for ii = 1:length(fields(stain_fname))
    for j = 1:length(stains)
        % Output filepath
        fpath = fullfile(dpath, stain_subid{ii}, 'stain');
        % Load the tissue mask
        mask = path_heatmap.(stain_subid{ii}).(stains{j}).mask;
        % Load the heatmap
        hm = path_heatmap.(stain_subid{ii}).(stains{j}).heatmap;
        % Figure Name
        stain_name = append(stain_fname.(stain_subid{ii}).(stains{j}),'_stain');
        figname = append(stain_name,'_heatmap_',num2str(cube_side));
        title_str = append(stain_subid{ii},' ', tstr{j});
        % Normalize each subvolume by largest dynamic range
        if strcmp(stains{j},'ab')
            hm = hm ./ ab_max;
            l = ab_min;
            u = ab_max;
        else
            hm = hm ./ pt_max;
            l = pt_min;
            u = pt_max;
        end
        % Plot heatmap
        plot_save_heatmap(1,hm,0,[0,1],mask,title_str,'(a.u.)',fpath,figname)
    end
end
%}

%% Vascular heatmap: length density, branch density, volume fraction
%%% Iterate over each subject
for ii = 1:length(subid)
    %% Generate Subvolumes
    %%% Load Data struct containing graph and angio (masked segmentation)
    graph = fullfile(dpath,subid{ii},subdir,segdir,graph_name);
    graph = load(graph);
    % Separate segmentation (angio) and the graph
    seg = graph.Data.angio;
    graph = graph.Data.Graph;   
    nodes = graph.nodes;

    %%% Load tissue mask
    mask = fullfile(dpath,subid{ii},mdir,'mask_tiss.mat');
    mask = load(mask);
    mask = mask.mask_tiss;
    
    %%% Create array of end node positions
    end_node_pos = nodes(graph.endNodes,:);
    % Calculate number of cubes in z dimension
    Nz = ceil(size(seg,3) ./ n_z);

    %% Iterate over the segmentation
    % Initialize the volume fraction heat map matrix 
    vf_mat = zeros(size(seg,1), size(seg,2), size(seg,3));
    ld_mat = zeros(size(seg,1), size(seg,2), size(seg,3));
    bd_mat = zeros(size(seg,1), size(seg,2), Nz);
    % Index to track the branch density matrix depth
    z_bd = 0;    
    % Iterate over the z-axis
    for z = 1:n_z:size(seg,3)
        z_bd = z_bd + 1;
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

                %% Calculate metrics if cube has segmentation
                % Check for segmentation within cube
                if sum(cube(:)) > 1
                    %% Volume Fraction
                    % Calculate volume fraction of cube
                    vf = sum(cube(:)) ./ cube_vol_vox;
                    % Assign to the vf matrix
                    vf_mat((x:xf), (y:yf), (z:zf)) = vf;
    
                    %% Branch Density
                    % Find branch nodes within the bounds of the cube
                    idx = find(...
                        x <= end_node_pos(:,2) &  end_node_pos(:,2) < xf &...
                        y <= end_node_pos(:,1) &  end_node_pos(:,1) < yf &...
                        z <= end_node_pos(:,3) &  end_node_pos(:,3) < zf);
                    % Extract number of branch points for each end node in cube
                    nb = sum(graph.nB(idx));
                    % Calculate branch density
                    bd = nb ./ cube_vol_um;
                    % Add branch density to the heatmap matrix
                    bd_mat((x:xf), (y:yf), z_bd) = bd;
                    
                    %% Length Density
                    %%% Create graph with subvolume (cube)
                    [g, ~] = seg_to_graph(cube, vox);
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
                        rm_loops_parallel(g.nodes, g.edges, cube,...
                                          delta, v_min, mv_iter, viz);
                    % Update graph with nodes and edges
                    g.nodes = nodes_rm;
                    g.edges = edges_rm;
                    
                    %%% If there is only 1 edge, then skip metadata. There
                    % is a bug in the graphing code that raises an error
                    % when there is only one edge. This is legacy code and
                    % difficult to debug. Instead, just calculate Euclidean
                    % distance of the lone edge.
                    if size(edges_rm,1) > 1
                        % Compute the graph metadata
                        [Data] = init_graph(g);
                        % Calculate length density with subvolume (cube)
                        seglen_um = Data.Graph.segInfo.segLen_um;
                        seglen_tot_um = sum(seglen_um);
                        ld = seglen_tot_um ./ cube_vol_um;
                    elseif size(edges_rm,1) == 1
                        % Euclidean Distance
                        n1 = nodes_rm(edges_rm(1),:) .* vox;
                        n2 = nodes_rm(edges_rm(2),:) .* vox;
                        d = sqrt( (n1(1) - n2(1)).^2 +...
                                  (n1(2) - n2(2)).^2 +...
                                  (n1(3) - n2(3)).^2);
                        ld = d ./ cube_vol_um;
                    else
                        % in this case there is no graph
                        ld = 0;
                    end
                    % Add length density to matrix
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
    heatmap.(sub).bd = bd_mat;
    sprintf('\nFinished subject %s\n',sub)
    
    %% Take maximum intensity projection of tissue mask  
    %%% Check to see if the mask matrix needs to be truncated
    if size(mask,1) > size(seg,1)
        mask = mask(1:size(seg,1),:,:);
    end
    if size(mask,2) > size(seg,2)
        mask = mask(:,1:size(seg,2),:);
    end

    % Initialize masks matrix
    masks = zeros(size(seg,1), size(seg,2), Nz);

    % Initialize zmin and zmax to index the tissue mask depths
    zmin = 1;
    zmax = zmin + n_z - 1;
    % Iterate over heatmap depths
    for z = 1:Nz
        % Take MIP
        masks(:,:,z) = max(mask(:,:,zmin:zmax),[],3);
        % Calculate z depth bounds for tissue mask
        zmin = zmin + n_z;
        zmax = min([(zmax + n_z), size(seg,3)]);
    end
    
    % Add mask to heatmap struct
    heatmap.(sub).mask = masks;

    %% Plot a figure for each frame of heatmap
    if viz_individual
        % Path to output heatmaps
        roi_dir = strcat('ROI_',num2str(cube_side));
        heatmap_dir = fullfile(mpath,'heatmaps',sub,roi_dir);
        if ~isfolder(heatmap_dir)
            mkdir(heatmap_dir);
        end
    
        % Iterate over depths in volume fraction heat map
        plot_save_heatmap(Nz, heatmap_vf, 0, [], masks,'Volume Fraction',...
            '(a.u.)', heatmap_dir,'heatmap_vf')
        % Iterate over depths in length density heat map
        plot_save_heatmap(Nz, heatmap_ld, 0, [], masks,'Length Density',...
            'Length (\mu) / Volume (\mu^3)',heatmap_dir,'heatmap_ld')
        % Iterate over depths in branch density heat map
        plot_save_heatmap(Nz, bd_mat, 0, [], masks,'Branch Density',...
            'Branches / Volume (\mu^3)',heatmap_dir,'heatmap_bd')
    end

    sprintf('FINISHED HEATMAP FOR SUBJECT %s',sub)
end

% Save the heatmap struct
heat_out = append('heatmap_',num2str(cube_side),'.mat');
heat_out = fullfile(mpath, heat_out);
save(heat_out,'heatmap','-v7.3');

%% Generate Heat Maps (not normalized)
% Individually generate a heat map for each subject. The colorbar scales
% are not normalied to one another.
%{
% Load Heat map
heatmap = load(fullfile(mpath, 'heatmap.mat'));
heatmap = heatmap.heatmap;
% Iterate over subject ID list
for ii = 1:length(subid)
    %%% Load the heatmap for the subject
    sub = subid{ii};
    heatmap_vf = heatmap.(sub).vf;
    heatmap_ld = heatmap.(sub).ld;
    heatmap_bd = heatmap.(sub).bd;
    masks = heatmap.(sub).mask;

    %%% Output filepath
    heatmap_dir = fullfile(mpath,'heatmaps',sub,roi_dir);
    if ~isfolder(heatmap_dir)
        mkdir(heatmap_dir);
    end

    % Iterate over depths in volume fraction heat map
    plot_save_heatmap(Nz, heatmap_vf, 0, [],masks,'Volume Fraction',...
        '(a.u.)', heatmap_dir,'heatmap_vf')
    % Iterate over depths in length density heat map
    plot_save_heatmap(Nz, heatmap_ld, 0, [],masks,'Length Density',...
        'Length (\mu) / Volume (\mu^3)',heatmap_dir,'heatmap_ld')
    % Iterate over depths in branch density heat map
    plot_save_heatmap(Nz, heatmap_bd, 0, [],masks,'Branch Density',...
        'Branches / Volume (\mu^3)',heatmap_dir,'heatmap_bd')
end
%}

%% Generate heat maps - normalized across subjects
% Iterate over each metric, choose one depth from each subject, normalize
% the colorbar across all subjects for this metric.

% Load Heat map
heat_out = append('heatmap_',num2str(cube_side),'.mat');
heat_out = fullfile(mpath, heat_out);
heatmap = load(heat_out);
heatmap = heatmap.heatmap;

% Vascular metrics (volume fraction, length density, branch density)
metrics = {'vf','ld','bd'};

%%% Identify minimum / maximum of each metric
% Set the maximum to be the 95th percentile of each metric. This will
% scale the colorbar to account for outliers.
vf = [];
ld = [];
bd = [];
for ii = 1:length(subid)
    vf = [vf, reshape(heatmap.(subid{ii}).vf(:),1,[])];
    ld = [ld, reshape(heatmap.(subid{ii}).ld(:),1,[])];
    bd = [bd, reshape(heatmap.(subid{ii}).bd(:),1,[])];
end
% Minimum of each metric
vf_min = min(vf);
ld_min = min(ld);
bd_min = min(bd);
% Maximum = 90th percentile of each metric
vf_max = prctile(vf,95);
ld_max = prctile(ld,95);
bd_max = prctile(bd,95);

%%% Generate normalized heatmaps
% Variable for whether or not to invert the heatmap
flip_cbar = 0;
% Iterate over each subject
for ii = 1:length(subid)
    %%% Load the heatmap for the subject
    sub = subid{ii};
    heatmap_vf = heatmap.(sub).vf;
    heatmap_ld = heatmap.(sub).ld;
    heatmap_bd = heatmap.(sub).bd;
    masks = heatmap.(sub).mask;

    %%% Output filepath for figures
    roi_dir = strcat('ROI_',num2str(cube_side));
    heatmap_dir = fullfile(mpath,'heatmaps',sub,roi_dir);
    if ~isfolder(heatmap_dir)
        mkdir(heatmap_dir);
    end

    %%% Plot first depth for each heatmap
    plot_save_heatmap([], heatmap_vf, flip_cbar, [vf_min, vf_max],...
        masks,'Volume Fraction','(a.u.)',...
        heatmap_dir,'rescaled_heatmap_vf')
    plot_save_heatmap([], heatmap_ld, flip_cbar, [ld_min, ld_max],...
        masks,'Length Density','Length (\mu) / Volume (\mu^3)',...
        heatmap_dir,'rescaled_heatmap_ld')
    plot_save_heatmap([], heatmap_bd, flip_cbar, [bd_min, bd_max],...
        masks,'Branch Density','Branches / Volume (\mu^3)',...
        heatmap_dir,'rescaled_heatmap_bd')
end

%% Plot and save the heat maps
function plot_save_heatmap(Ndepths, heatmaps, flip_cbar, colorbar_range,...
    masks, tstr, cbar_label, dpath, fname)
% PLOT_SAVE_HEATMAP use imagesc and set background = 0
% INPUT
%   Ndepths (int): number of depths in z dimension
%   heatmaps (double matrix): heatmaps of vascular metric
%   flip_cbar (logical): reverse the direction of the colorbar
%   colorbar_range (double array): [min, max]
%   masks (double): tissue mask (1=tissue, 0=other)
%   tstr (string): figure title
%   cbar_label (string): colorbar label
%   dpath (string): data directory path
%   fname (string): name of figure to save

%%% Set the number of depths to iterate for each heatmap
% If the number of depths is not specified, then set it equal to the number
% of z dimensions.
if isempty(Ndepths)
    Ndepths = size(heatmaps,3);
end
% Set fontsize for the heatmap figure
fontsize = 40;

%%% Iterate over frames in z dimension
for d = 1:Ndepths
    %%% Heatmap of j_th frame from the length density
    fh = figure();
    fh.WindowState = 'maximized';
    % If there are multiple heatmaps in the matrix
    if size(heatmaps,3) > 1
        heatmap = heatmaps(:,:,d);
    % Here it is just a single frame of a heatmap
    else
        heatmap = heatmaps;
    end
    % Initialize heatmap
    h = imagesc(heatmap);

    %%% Initialize colorbar
    % If the colorbar_range is passed in, then extract min & max
    if ~isempty(colorbar_range)
        cmap_min = colorbar_range(1);
        cmap_max = colorbar_range(2);
    % Otherwise, set limits from the current heatmap
    else
        % Min = Lowest value also greater than zero
        cmap_min = min(heatmap(heatmap(:)>0));
        % Max = Find the 95th percentile for upper limit
        cmap_max = prctile(heatmap(:), 95);
    end
    % Initialize the colormap limits
    cmap = jet(256);
    clim(gca, [cmap_min, cmap_max]);
    % Initialize colormap and colorbar
    if flip_cbar
        colormap(flipud(cmap));
    else
        colormap(cmap);
    end
    c = colorbar;

    %%% Apply tissue mask to the heatmap to remove background
    alpha_mask = double(masks(:,:,d));
    set(h, 'AlphaData', alpha_mask);

    %%% Configure figure parameters
    % Update title string with depth
    if size(heatmaps,3) > 1
        title_str = append(tstr, ' Depth ', num2str(d));
    else
        title_str = tstr;
    end
    title(title_str,'Interpreter','none');
    set(gca, 'FontSize', fontsize);
    % Label the colorbar    
    c.Label.String = cbar_label;
    % Offset colorbar label to the right of colorbar
    c.Label.Position = [10 (cmap_max - (cmap_max-cmap_min)/2)];
    c.Label.Rotation = 270;
    % Increase fontsize of colorbar
    c.FontSize = 40;
    % Remove x and y tick labels
    set(gca,'Yticklabel',[]);
    set(gca,'Xticklabel',[]);
    
    %%% Save figure as PNG
    % If there are multiple heatmaps in the matrix, save vascular heatmap
    if size(heatmaps,3) > 1
        fout = append(fname, '_', num2str(d));
    % Otherwise, save the pathology heatmap
    else
        fout = fname;
    end
    % If the colorbar is reversed then add suffix to filename
    if flip_cbar
        fout = append(fout, '_flip_cbar');
    end

    fout = fullfile(dpath, fout);
    saveas(gca, fout,'png');
    close;

end
end

%% Pathology heatmap
function [hm] = pathology_heatmap(pix, square_side, stain, th, mask)
%PATHOLOGY_HEATMAP generate heatmap with pathology
% Divide the pathology stain into isotropic squares, calculate the density
% of each square, add the density to a matrix.
%
% INPUTS:
%   pix (double array): pixel size (x,y) (microns)
%   square_side (uint): heatmap square dimension (microns)
%   stain (uint8 matrix): pathology stain
%   th (double): minimum threshold for segmenting AB or AT8 (after taking
%               complement of the staining)
%   mask (logical matrix): tissue mask for the pathology stain
% OUTPUTS:
%   hm (double matrix): masked heatmap of the pathology stain

%%% Verify that the stained image is scaled between [0,1]
% The inversion assumes that the image is scaled between [0,1], so this
% will return the incorrect value otherwise.
r = range(stain);
assert(0<=min(r), 'The stained image is not scaled between [0,1]');
assert(max(r)<=1, 'The stained image is not scaled between [0,1]');

%%% Translate cube size to pixel size
% Number of pixels in each dimension to create square
nx = floor(square_side ./ pix(1));
ny = floor(square_side ./ pix(2));

%%% Complement image and segment plaques (minimum threshold)
% Complement the image (you look nice, today)
stain = imcomplement(stain);
% Set pixels above background intensity equal to 1
stain(stain > th) = 1;
% Mask the stain to remove non-tissue pixels
stain(~mask) = 0;
% Set pixels < 1 equal to zero
stain(stain ~= 1) = 0;

%%% Generate heat map
% Initialize heatmap matrices for A-beta and p-tau
hm = zeros(size(stain,1), size(stain,2));
% Iterate over rows
for x = 1:nx:size(stain,1)
    % Iterate over columns
    for y = 1:ny:size(stain,2)
        %%% Crop segmentation into isotropic square
        % Initialize end indices for each axis
        xf = x + nx - 1;
        yf = y + ny - 1;
        % Take minimum of matrix dimensions and end indices
        xf = min(xf, size(stain,1));
        yf = min(yf, size(stain,2));
        % Take cube from segmentation
        path_square = stain((x:xf), (y:yf));

        %%% Calculate average plaque intensity in subvolume
        hm((x:xf), (y:yf)) = mean(path_square(:));
    end
end

% Mask the pathology heatmap with tissue mask
hm(~mask) = 0;
end

%% Read TIFF and extract image properties
function [stain, mask, res] = import_pathology(dpath, stain_name, mask_name)
%IMPORT_PATHOLOGY read the TIF file for the pathology (A-beta or p-tau)
% Retrieve TIFF metadata for the staining & compute the resolution of each
% pixel. Then, load the TIFF and the respective tissue mask. Rescale the
% staining to [0,1] to standardize across all images. Then, apply the
% tissue mask.
% INPUTS:
%   dpath (string): directory path to the staining TIFF file
%   stain_name (string): filename of staining TIFF file
%   mask_name (string): filename of tissue mask TIFF file
% OUTPUTS:
%   stain (double matrix): masked staining, scaled between [0,1]
%   mask (logical matrix): tissue mask for the pathology slide
%   bg (double matrix): background patch for normalizing the pathology.
%                       This is a region of the stain from the white matter
%                       that does not include any plaques.
%   res (double array): resolution of a pixel [x,y]

%%% Disable Tiff read warnings for unknown tags
id = 'imageio:tiffmexutils:libtiffWarning';
warning('off',id);

%%% Retrieve metadata and calculate resolution
% Read staining TIFF info
stain_file = append(stain_name, '.tif');
stain_file = fullfile(dpath, stain_file);
stain_metadata = imfinfo(stain_file);
% X and Y resolution (pixels / micron)
xres = stain_metadata.XResolution;
yres = stain_metadata.YResolution;
% Invert resolution (microns / pixel)
xres = 1./xres;
yres = 1./yres;
res = [xres, yres];

%%% Load TIFF file
% Read staining TIFF (RGB) & convert to grayscale
stain = Tiff(stain_file,'r');
stain = read(stain);
stain = rescale(stain, 0, 1);
stain = rgb2gray(stain);
% Load tissue mask
mask_fname = append(mask_name, '.tif');
mask_file = fullfile(dpath, mask_fname);
mask = Tiff(mask_file,'r');
mask = logical(read(mask));
end
