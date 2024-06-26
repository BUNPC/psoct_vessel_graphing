%% Apply tissue masks to the heatmaps and generate distributions
% This should be run after the script "generate_heatmap_subgraphs.m"
% This script will save "heatmap_[#]_distro.mat" where # corresponds to
% the ROI size. This will be saved to the metrics path.
% The next script in the sequence is
% "vessel_analyses_wm_gm_sulci_gyri_main.m"

% TODO:
% - load the heatmaps + masks
% - For each subject:
%   - remove ROIs outside of tissue boundaries
%   - save the ROIs within tissue mask to the struct hm_distro
% - Output the struct "heatmap_[#]_distro.mat"

%% Clear workspace & add top-level directory
clear; clc; close all;

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

%% Define directory paths
% Path to top-level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
   
% Volume and graph directories
voldir = '/dist_corrected/volume/ref/';
maskdir = 'dist_corrected/volume/ref/masks/';
heatmap_dir = 'metrics/gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/';

% Metrics output path
mpath = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
    'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];

% Masks corresponding to each respective graph
masks = {'mask_tiss.mat','mask_gyri','mask_sulci','mask_gm','mask_wm'};

% IDs of each subject
subid = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_7126',...
         'NC_6839',  'NC_6974',  'NC_8653',  'NC_21499', 'NC_8095'};

%%% Heatmap struct to import
% Isotropic cube length (microns)
cube_side = 1000;
% create name of matlab struct
hm = append('heatmap_',num2str(cube_side),'.mat');
% Voxel size (microns)
vox = [12, 12, 15];
% Compute number of voxels in x,y dimensions for each cube
n_x = floor(cube_side ./ vox(1));
n_y = floor(cube_side ./ vox(2));
n_z = floor(cube_side ./ vox(3));

% struct for storing heatmap distributions for each subject
hm_distro = struct();

% Load the heatmap struct
hm = load(fullfile(dpath,heatmap_dir,hm));
hm = hm.heatmap;

%% Iterate through subjects
for ii = 1:length(subid)
    % Specify the subject ID
    sub = subid{ii};
    fprintf('Calculating Metrics for Subject %s\n', sub);
    
    % Load vascular metrics heatmaps
    vf = hm.(sub).vf;
    ld = hm.(sub).ld;
    bd = hm.(sub).bd;

    % Create struct for storing masks
    masks = struct();

    % Load all premade masks (tissue, gyri, sulci, GM, WM)
    mask_tiss = load(fullfile(dpath,sub,maskdir,'mask_tiss.mat'),'mask_tiss');
    mask_gyri = load(fullfile(dpath,sub,maskdir,'mask_gyri.mat'));
    mask_sulci = load(fullfile(dpath,sub,maskdir,'mask_sulci.mat'));
    mask_gm = load(fullfile(dpath,sub,maskdir,'mask_gm.mat'));
    mask_wm = load(fullfile(dpath,sub,maskdir,'mask_wm.mat'));
    mask_tiss = mask_tiss.mask_tiss;
    mask_gyri = mask_gyri.mask_gyri;
    mask_sulci = mask_sulci.mask_sulci;
    mask_gm = mask_gm.mask_gm;
    mask_wm = mask_wm.mask_wm;

    %%% Find the maximum depth to use for the masks
    zmin = min([size(mask_tiss,3), size(mask_gyri,2), size(mask_sulci,3),...
        size(mask_gm,3), size(mask_wm,3)]);

    % Truncate each volume to the maximum depths
    masks.tiss =     mask_tiss(:, :, 1:zmin);
    masks.wm =       mask_wm(:, :, 1:zmin);
    masks.gm =       mask_gm(:, :, 1:zmin);
    masks.sulci =    mask_sulci(:, :, 1:zmin);
    masks.gyri =     mask_gyri(:, :, 1:zmin);

    % Create remaining masks (WM & GM partitions of sulci/gyri)
    masks.gm_sulci = masks.gm .* masks.sulci;
    masks.wm_sulci = masks.wm .* masks.sulci;
    masks.gm_gyri = masks.gm .* masks.gyri;
    masks.wm_gyri = masks.wm .* masks.gyri;

    %% Perform MIP of each mask at specified ROI z dimension
    % The tissue masks are for the entire depth, but the heatmaps are
    % already averaged at a predetermined depth for each ROI.
    
    % struct to store the MIP masks
    masks_mip = struct();

    % Number of cubes in z dimension from one of the heatmaps
    Nz = size(vf,3);

    % Iterate over the masks
    f = fields(masks);
    for j=1:length(f)
        % Load the respective mask for this iteration
        mask = masks.(f{j});
        % Initialize matrix to store maximum intensity projection of mask
        mask_mip = zeros(size(mask,1), size(mask,2), Nz);
        
        %%% MIP at the specific depth
        % Initialize zmin and zmax to index the tissue mask depths
        zmin = 1;
        zmax = zmin + n_z - 1;
        % Iterate over heatmap depths
        for z = 1:Nz
            % Take MIP
            mask_mip(:,:,z) = max(mask(:,:,zmin:zmax),[],3);
            % Calculate z depth bounds for tissue mask
            zmin = zmin + n_z;
            zmax = min([(zmax + n_z), size(mask,3)]);
        end
        % Save the mask
        masks_mip.(f{j}) = mask_mip;
    end

    %% Iterate through masks for each tissue region
    % Apply the mask
    % Check if the cubic grid contains the mask
    %   - if so, then store the value from this cube

    f = fields(masks);
    for j=1:length(f)
        %%% Array to store single value from each cube within mask
        % Calculate maximum number of cubes in tissue volume. The z
        % dimension was already compressed in the prior step when creating
        % the heatmap
        Nz = size(vf,3);
        Ny = ceil(size(vf,1) ./ n_y);
        Nx = ceil(size(vf,2) ./ n_x);
        n_cubes = Nz * Ny * Nx;
        % Initialize array for each metric
        vf_distro = zeros(n_cubes,1);
        ld_distro = zeros(n_cubes,1);
        bd_distro = zeros(n_cubes,1);

        %%% Load the respective mask for this iteration
        mask = masks_mip.(f{j});
        % Determine whether the mask matrix needs to be truncated
        if size(mask,1) > size(vf,1)
            mask = mask(1:size(vf,1),:,:);
        end
        if size(mask,2) > size(vf,2)
            mask = mask(:,1:size(vf,2),:);
        end
        % Retain the non-zero voxels within the mask
        idx = find(mask);
        % Apply mask to each heatmap
        vf_masked = vf(idx) .* double(mask(idx));
        ld_masked = ld(idx) .* double(mask(idx));
        bd_masked = bd(idx) .* double(mask(idx));

        %% Iterate through cubic grid and retain values within mask
        % Counter for storing the distribution values
        cnt = 1;
        % Iterate over z dimension (already averaged from last step)
        for z = 1:size(vf,3)
            % Iterate over rows
            for x = 1:n_x:size(vf,1)
                % Iterate over columns
                for y = 1:n_y:size(vf,2)
                %% Crop segmentation into cube
                % Initialize end indices for each axis
                xf = x + n_x - 1;
                yf = y + n_y - 1;
                % Take minimum of matrix dimensions and end indices
                xf = min(xf, size(vf,1));
                yf = min(yf, size(vf,2));
                % Take cube from mask
                m = mask((x:xf), (y:yf), z);
                % Determine if the cube contains any portion of mask
                if any(m)
                    % Retain single value from cube
                    vf_distro(cnt) = max(vf((x:xf), (y:yf), z),[],"all");
                    ld_distro(cnt) = max(ld((x:xf), (y:yf), z),[],"all");
                    bd_distro(cnt) = max(bd((x:xf), (y:yf), z),[],"all");
                    % Iterate the counter
                    cnt = cnt + 1;
                end
                end
            end
        end
        % Store distribution
        hm_distro.(sub).(f{j}).vf = vf_distro(1:cnt);
        hm_distro.(sub).(f{j}).ld = ld_distro(1:cnt);
        hm_distro.(sub).(f{j}).bd = bd_distro(1:cnt);
    end

    fprintf('Finished Subject %s\n', sub);
end

%% Save output
% Add the ROI size to the output filename
fout = append('heatmap_distro_',num2str(cube_side),'.mat');
save(fullfile(mpath, fout), 'hm_distro','-v7.3');