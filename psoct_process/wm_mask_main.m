%% Create a mask for the white matter
% Use the scattering map to segment just the white matter

clear; clc; close all;

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

%% Initialize subject ID lists
%%% All subjects to analyze
subid = {'AD_10382', 'AD_20832', 'AD_20969',...
         'AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912',...
         'CTE_7019', 'CTE_8572','CTE_7126',...
         'NC_6839',  'NC_6974', 'NC_8653',...
         'NC_21499', 'NC_301181'};
subid = {'AD_21354'};

%%% stacks that require truncating
trunc = {'AD_10382','AD_20832','CTE_7126','CTE_8572',...
         'NC_6839',  'NC_6974', 'NC_8653',...
         'NC_21499', 'NC_301181'};
% Array of final image in stack
zmins = [199, 229, 201, 198, 181, 198, 187, 178, 220];
% Create dictionary to store last image in stack
d = dictionary(trunc, zmins);

%% Initialize directories and filenames

%%% Directories 
% Upper level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T';    
% Subfolder with normalized volume
subdir = '/dist_corrected/volume';
% Subfolder containing ref files
subdir1 = '/dist_corrected/volume/ref';
% Combined segmentations subfolder
subdir2 = 'combined_segs';
% Combined segmentation subfolder
segdir = 'gsigma_2-3-4_3-5-7_5-7-9_7-9-11';
% Scattering subdirectory
scatdir = 'fitting_4x/';

%%% Filenames
% Combined segmentations (.MAT)
seg_name = 'seg';
% non-normalized PSOCT volume
vol_name = 'ref';
% Normalized PSOCT tissue volume
voln_name = 'ref_4ds_norm_inv';

%% Iterate through subjects. Create and apply mask.
for ii = 1:length(subid)
    %%% Debugging information
    fprintf('---------Starting Subject %s---------\n',subid{ii})

    %%% Import volumes, mask, and segmentation
    % Import normalized volume
    fprintf('importing normalized volume\n')
    fpath = fullfile(dpath, subid{ii}, subdir, strcat(voln_name,'.tif'));
    voln = TIFF2MAT(fpath);   
    
    % Import non-normalized volume
    fprintf('importing non-normalized volume\n')
    fpath = fullfile(dpath, subid{ii}, subdir1, strcat(vol_name,'.mat'));
    vol = load(fpath); 
    vol = vol.vol;
    
    % Import combined segmentation file
    fprintf('importing segmentation\n')
    fpath = fullfile(dpath, subid{ii}, subdir, subdir2, segdir, 'seg.mat');
    seg = load(fpath,'seg');
    seg = seg.seg;

    % Import scattering Map
    fprintf('importing mask\n')
    fpath = fullfile(dpath, subid{ii}, scatdir, 'mus1.tif');
    mus = TIFF2MAT(fpath);

    %% Register scattering map to first slice
    slice = vol(:,:, 1);
    [optimizer, metric]  = imregconfig('monomodal');
    mus_reg = imregister(mus,slice,'translation',optimizer,metric);
    figure;
    imshowpair(slice,mus_reg)

    %% Threshold the mus
    % Tissue foreground mask
    fg = mus_reg;
    fg(fg>0) = 1;
    se = strel('disk', 10);
    fg = imerode(fg, se);
    figure; imshow(fg);
    
    mus_th = mus_reg;
    mus_th(mus_th<=15) = 0;
    mus_th = bwareaopen(mus_th, 50);
    figure; imagesc(mus_th);
    
    % Close the holes
    se = strel('disk', 10);
    close_mus = imclose(mus_th,se);
    figure; imagesc(close_mus);

    % Remove voxels with fewer than 20 connections
    wm_mask = bwareaopen(close_mus, 20);
    figure; imagesc(wm_mask);

    % white matter mask
    figure; imshow(wm_mask);

    % Gray matter mask & foreground mask
    gm_mask = ~wm_mask .* fg;
    figure; imshow(gm_mask);

    %% Apply WM mask to tissue
    wm = slice .* wm_mask;
    gm = slice .* gm_mask;

    figure; imshow(wm); title('white matter')
    figure; imshow(gm); title('gray matter')

    %% Set zmin if truncating finals slices
    % Truncate the mask, and the others will be subsequently truncated.
    sid = subid{ii};
    if isKey(d, {sid})
        % Retrieve the zmin from the array
        zmin = d({sid});
        % Truncate the mask
        mask = mask(:,:,1:zmin);
    end

    %%% Verify dimensions of mask and volumes
    % The segmentation file was made using the normalized volume. Therefore
    % the dimensions of the mask must be made equal to those of the seg.
    if any(size(mask) ~= size(vol)) || any(size(mask) ~= size(voln))
        %%% Truncate mask & vol (ref) to match dimensions of normalized
        m = size(mask);
        vn = size(voln);
        % Find smallest dimensions
        ymin = min(m(1), vn(1));
        xmin = min(m(2), vn(2));
        zmin = min(m(3), vn(3));
        % Truncate all volumes
        mask = mask(1:ymin, 1:xmin, 1:zmin);
        vol = vol(1:ymin, 1:xmin, 1:zmin);
        voln = voln(1:ymin, 1:xmin, 1:zmin);
        seg = seg(1:ymin, 1:xmin, 1:zmin);
        % Assert that the new matrices are all equivalent
        assert(isequal(size(mask), size(voln)),'Mask and voln dimensions mismatch')
        assert(isequal(size(mask), size(vol)),'Mask and voln dimensions mismatch')
    end

    %%% Apply mask
    fprintf('apply mask\n')
    try
        % Mask the non-normalized volume (volm)
        volm = vol .* mask;
        % Mask the normalized volume (volnm)
        volnm = voln .* uint16(mask);
        % Mask the segmentation (segm)
        segm = seg .* uint8(mask);
        % Clear variables to save memory
        clear vol; clear voln; clear seg;    
    catch
        %%% catch where mask if different data type
        fprintf('failed to apply mask on subject %s\n',sid)
        continue
    end
    
    %%% Save masked volumes
    fprintf('saving masked volume\n')
    % Export masked non-normalized volume (volm)
    volm_out = fullfile(dpath, subid{ii}, subdir1,...
        strcat(vol_name,'_refined_masked.tif'));
    segmat2tif(volm, volm_out);
    
    % Export masked normalized volume (volnm)
    volnm_out = fullfile(dpath, subid{ii}, subdir1,...
        strcat(voln_name,'_refined_masked.tif'));
    segmat2tif(volnm, volnm_out);

    %%% Debugging Info
    fprintf('---------Finished Subject %s---------\n',subid{ii})
end
