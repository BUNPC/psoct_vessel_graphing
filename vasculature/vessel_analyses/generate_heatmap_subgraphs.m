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
% Set maximum number of threads equal to number of threads for script
NSLOTS = str2num(getenv('NSLOTS'));
maxNumCompThreads(NSLOTS);

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
% Graph (without loops) from masked segmentation
graph_name = 'seg_refined_masked_rmloop_graph_data.mat';

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

%%% Struct for storing heat map
heatmap = struct();

%% Generate grid - average length density, branch density, volume fraction

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
    vf_mat = zeros(size(seg,1), size(seg,2), size(seg,3));
    ld_mat = zeros(size(seg,1), size(seg,2), size(seg,3));
    bd_mat = zeros(size(seg,1), size(seg,2), Nz);
    % Index to track the branch density matrix depth
    z_bd = 0;
    % Calculate the size of each cube (in voxels)
    cube_vol_vox = n_x * n_y * n_z;
    cube_vol_um = cube_vol_vox * vox(1) * vox(2) * vox(3);
    
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

                %% Calculate length density if cube has segmentation
                % Check for segmentation within cube
                if sum(cube(:)) > 1
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
    % Path to output heatmaps
    heatmap_out = fullfile(mpath,'heatmaps',sub);
    if ~isfolder(heatmap_out)
        mkdir(heatmap_out);
    end

    % Iterate over depths in volume fraction heat map
    plot_save_heatmap(Nz, heatmap_vf, masks,'Volume Fraction',...
        '(a.u.)', heatmap_out,'heatmap_vf')
    % Iterate over depths in length density heat map
    plot_save_heatmap(Nz, heatmap_ld, masks,'Length Density',...
        'Length (\mu) / Volume (\mu^3)',heatmap_out,'heatmap_ld')
    % Iterate over depths in branch density heat map
    plot_save_heatmap(Nz, bd_mat, masks,'Branch Density',...
        'Branches / Volume (\mu^3)',heatmap_out,'heatmap_bd')
end


%% Save the heatmap
heat_out = fullfile(mpath, 'heatmap.mat');
save(heat_out,'heatmap','-v7.3');

%% Generate Heat Maps
% Individually generate a heat map for each subject. The colorbar scales
% are not normalied to one another.

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
    heatmap_out = fullfile(mpath,'heatmaps',sub);
    if ~isfolder(heatmap_out)
        mkdir(heatmap_out);
    end

    % Iterate over depths in volume fraction heat map
    plot_save_heatmap(Nz, heatmap_vf,[],masks,'Volume Fraction',...
        '(a.u.)', heatmap_out,'heatmap_vf')
    % Iterate over depths in length density heat map
    plot_save_heatmap(Nz, heatmap_ld,[],masks,'Length Density',...
        'Length (\mu) / Volume (\mu^3)',heatmap_out,'heatmap_ld')
    % Iterate over depths in branch density heat map
    plot_save_heatmap(Nz, heatmap_bd,[],masks,'Branch Density',...
        'Branches / Volume (\mu^3)',heatmap_out,'heatmap_bd')
end

%% Genearte heat maps - normalized across subjects
% Iterate over each metric, choose one depth from each subject, normalize
% the colorbar across all subjects for this metric.

% Load Heat map
heatmap = load(fullfile(mpath, 'heatmap.mat'));
heatmap = heatmap.heatmap;

% Vascular metrics (volume fraction, length density, branch density)
metrics = {'vf','ld','bd'};

%%% Identify minimum / maximum of each metric
% Initialize the minimum values
vf_min = min(heatmap.(subid{1}).vf(:));
ld_min = min(heatmap.(subid{1}).ld(:));
bd_min = min(heatmap.(subid{1}).bd(:));
% Initialize the maximum values
vf_max = max(heatmap.(subid{1}).vf(:));
ld_max = max(heatmap.(subid{1}).ld(:));
bd_max = max(heatmap.(subid{1}).bd(:));
% Iterate over subject ID list
for ii=1:length(subid)
    % Volume fraction
    vf_min = min(vf_min, min(heatmap.(subid{ii}).vf(:)));
    vf_max = max(vf_max, max(heatmap.(subid{ii}).vf(:)));
    % Length Density
    ld_min = min(ld_min, min(heatmap.(subid{ii}).ld(:)));
    ld_max = max(ld_max, max(heatmap.(subid{ii}).ld(:)));
    % Branch Density
    bd_min = min(bd_min, min(heatmap.(subid{ii}).bd(:)));
    bd_max = max(bd_max, max(heatmap.(subid{ii}).bd(:)));
end

%%% Generate normalized figures
% Set number of depths equal to 1
Nz = 1;

%% Iterate over subject ID list
for ii = 1:length(subid)
    %%% Load the heatmap for the subject
    sub = subid{ii};
    heatmap_vf = heatmap.(sub).vf;
    heatmap_ld = heatmap.(sub).ld;
    heatmap_bd = heatmap.(sub).bd;
    masks = heatmap.(sub).mask;

    %%% Output filepath for figures
    heatmap_out = fullfile(mpath,'heatmaps',sub);
    if ~isfolder(heatmap_out)
        mkdir(heatmap_out);
    end

    %%% Plot first depth for each heatmap
    plot_save_heatmap(1, heatmap_vf, [vf_min, vf_max], masks,...
        'Volume Fraction','(a.u.)',...
        heatmap_out,'rescaled_heatmap_vf')
    plot_save_heatmap(1, heatmap_ld, [ld_min, ld_max], masks,...
        'Length Density','Length (\mu) / Volume (\mu^3)',...
        heatmap_out,'rescaled_heatmap_ld')
    plot_save_heatmap(1, heatmap_bd, [bd_min, bd_max], masks,...
        'Branch Density','Branches / Volume (\mu^3)',...
        heatmap_out,'rescaled_heatmap_bd')
end




%% Plot and save the heat maps
function plot_save_heatmap(Ndepths, heatmaps, colorbar_range, masks,...
    tstr, ystr, dpath, fname)
% PLOT_SAVE_HEATMAP use imagesc and set background = 0
% INPUT
%   Ndepths (int): number of depths in z dimension
%   heatmaps (double matrix): heatmaps of vascular metric
%   colorbar_range (double array): [min, max]
%   masks (double): tissue mask (1=tissue, 0=other)
%   tstr (string): figure title
%   ystr (string): colorbar label
%   dpath (string): data directory path
%   fname (string): name of figure to save

fontsize = 40;

%%% Iterate over frames in z dimension
for d = 1:Ndepths
    %%% Heatmap of j_th frame from the length density
    fh = figure();
    fh.WindowState = 'maximized';
    heatmap = heatmaps(:,:,d);
    h = imagesc(heatmap);

    %%% Initialize colorbar
    % If the colorbar_range is passed in, then extract min & max
    if exist('colorbar_range', 'var')
        cmap_min = colorbar_range(1);
        cmap_max = colorbar_range(2);
    % Otherwise, set limits from the current heatmap
    else
        cmap_min = min(heatmap(heatmap(:)>0));
        cmap_max = max(heatmap(:));
    end
    % Initialize the colormap limits
    cmap = jet(256);
    clim(gca, [cmap_min, cmap_max]);
    % Initialize colormap and colorbar
    colormap(cmap);
    c = colorbar;

    %%% Find values equal to zero in heatmap
    % Take inverse to mask out the zero pixels
    alpha_mask = double(masks(:,:,d));
    set(h, 'AlphaData', alpha_mask);

    %%% Configure figure parameters
    % Update title string with depth
    title_str = append(tstr, ' Depth ', num2str(d));
    title(title_str);
    set(gca, 'FontSize', fontsize);
    % Label the colorbar
    yh = ylabel(c, ystr,"FontSize",fontsize,'Rotation',270);
    % Offset colorbar label to the right
    yh.Position(1) = yh.Position(1) + 0.2;
    % Increase fontsize of colorbar
    c.FontSize = 40;
    % Remove x and y tick labels
    set(gca,'Yticklabel',[]);
    set(gca,'Xticklabel',[]);
    
    %%% Save figure as PNG
    fout = append(fname, '_', num2str(d));
    fout = fullfile(dpath, fout);
    saveas(gca, fout,'png');
    close;

end
end