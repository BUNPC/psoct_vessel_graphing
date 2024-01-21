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
subid = {'AD_10382'};

%%% stacks that require truncating
trunc = {'AD_10382','AD_20832','CTE_7126','CTE_8572',...
         'NC_6839',  'NC_6974', 'NC_8653',...
         'NC_21499', 'NC_301181'};
% Array of final image in stack
zmins = [199,229, 201, 198, 181, 198, 187, 178, 220];
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
    
    % Import mask matrix
    fprintf('importing mask\n')
    fpath = fullfile(dpath, subid{ii}, subdir1, 'maskec.mat');
    mask = load(fpath);
    mask = mask.mask;
    % Convert mask to logical for matrix operations
    mask = imbinarize(mask);

    % Import combined segmentation file
    fprintf('importing segmentation\n')
    fpath = fullfile(dpath, subid{ii}, subdir, subdir2, segdir, 'seg.mat');
    seg = load(fpath,'seg');
    seg = seg.seg;
    
    %%% Set zmin if truncating finals slices
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
    
    %%% Overlays the segmentation and normalized masked volume
    fprintf('overlaying segmentation and masked normalized volume\n')
    fout = fullfile(dpath, subid{ii}, subdir, subdir2, segdir,...
                    strcat(seg_name, '_refined_mask_overlay_norm.tif'));
    overlay_vol_seg(volnm, segm, 'green', fout);
    
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
    
    % Export masked segmentation (segm)
    segm_out = fullfile(dpath, subid{ii}, subdir, subdir2, segdir,...
        strcat(seg_name,'_refined_masked.mat'));
    save(segm_out,'segm','-v7.3')
    segm_out = fullfile(dpath, subid{ii}, subdir, subdir2, segdir,...
        strcat(seg_name,'_refined_masked.tif'));
    segmat2tif(segm, segm_out);
    
    %%% Convert to graph
    fprintf('Generating graph data for sub %s\n',subid{ii})
    vox_dim = [12, 12, 15];
    fullpath = fullfile(dpath, subid{ii}, subdir, subdir2, segdir);
    fname_seg = strcat(seg_name,'_refined_masked');
    viz = false;
    rmloop_bool = false;
    seg_graph_init(segm, vox_dim, fullpath, fname_seg, viz, rmloop_bool)

    %%% Debugging Info
    fprintf('---------Finished Subject %s---------\n',subid{ii})
end
