%% Apply the tissue, white matter, and gray matter masks to segmentation
% This script applies the tissue mask to the segmentation. Then, it applies
% the white matter mask and then the inverse of the white matter mask.
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
         'CTE_6489', 'CTE_6912','CTE_7019','CTE_7126','CTE_8572',...
         'NC_6839',  'NC_6974', 'NC_8653','NC_21499', 'NC_301181'};

%%% stacks that require truncating
% Last depth to retain for each stack
zmins = [187, 165, 242, 165, 220,...
        242, 110, 198, 198, 198,...
        176, 198, 165, 165, 220];
% Create dictionary to store last image in stack
d = dictionary(subid, zmins);

%% Initialize directories and filenames

%%% Directories 
% Upper level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T';    
% Subfolder with normalized volume
subdir = '/dist_corrected/volume';
% Subfolder containing ref files
subdir1 = '/dist_corrected/volume/ref';
% Combined segmentation subfolder
segdir = '/combined_segs/gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/';
% Mask subfolder
mdir = '/dist_corrected/volume/ref/masks';

%%% Filenames
% Combined segmentations (.MAT)
seg_name = 'seg';
% non-normalized PSOCT volume
vol_name = 'ref';
% Normalized PSOCT tissue volume
voln_name = 'ref_4ds_norm_inv';

%% Iterate through subjects. Create and apply mask.
% parfor (ii = 2:length(subid),NSLOTS)
for ii = 15:length(subid)
    %%% Debugging information
    fprintf('\n---------Starting Subject %s---------\n',subid{ii})

    %% Import volumes, mask, and segmentation
    %%% Import normalized volume
    fpath = fullfile(dpath, subid{ii}, subdir, strcat(voln_name,'.tif'));
    voln = TIFF2MAT(fpath);   

    %%% Import non-normalized volume
    fpath = fullfile(dpath, subid{ii}, subdir1, strcat(vol_name,'.mat'));
    vol = load(fpath); 
    vol = vol.vol;

    %%% Import combined segmentation file
    fpath = fullfile(dpath, subid{ii}, subdir, segdir, 'seg.mat');
    seg = load(fpath,'seg');
    field = fieldnames(seg);
    field = field{:};
    seg = seg.(field);

    %%% Import tissue mask matrix
    fpath = fullfile(dpath, subid{ii}, subdir1, 'maskec.mat');
    mask = load(fpath);
    mask = mask.mask;
    % Convert mask to logical for matrix operations
    if ~isa(mask,'logical')
        mask = imbinarize(mask);
    end

    %%% Import white matter mask
    % Import the cleaned mask, if it exists
    wmask_path = fullfile(dpath, subid{ii}, mdir,'wm_mask.mat');
    wmaskec_path = fullfile(dpath, subid{ii}, mdir,'wm_maskec.mat');
    if isfile(fullfile(wmaskec_path))
        wm_mask = load(wmaskec_path);
    else
        wm_mask = load(wmask_path);
    end
    % Convert mask to logical for matrix operations
    field = fieldnames(wm_mask);
    field = field{:};
    wm_mask = wm_mask.(field);
    % Convert to logical if not already
    if ~isa(wm_mask,'logical')
        wm_mask = imbinarize(wm_mask);
    end
    
    %% Verify dimensions of masks and volumes
    % The segmentation file was made using the normalized volume. Therefore
    % the dimensions of the mask must be made equal to those of the seg.

    %%% Set zmin for truncating finals slices
    % Truncate the mask, and the others will be subsequently truncated.
    sid = subid{ii};
    % Retrieve the zmin from the array
    zmin = d({sid});
    % Truncate the tissue mask
    mask = mask(:,:,1:zmin);

    %%% Truncate mask & vol (ref) to match dimensions of normalized
    % Find smallest dimensions
    ymin = min([size(mask,1), size(wm_mask,1), size(voln,1)]);
    xmin = min([size(mask,2), size(wm_mask,2), size(voln,2)]);
    zmin = min([size(mask,3), size(wm_mask,3), size(voln,3)]);
    % Truncate all volumes
    mask =      mask(1:ymin, 1:xmin, 1:zmin);
    wm_mask =   wm_mask(1:ymin, 1:xmin, 1:zmin);
    vol =       vol( 1:ymin, 1:xmin, 1:zmin);
    voln =      voln(1:ymin, 1:xmin, 1:zmin);
    seg =       seg( 1:ymin, 1:xmin, 1:zmin);
    % Assert that the new matrices are all equivalent
    assert(isequal(size(mask),size(voln)),...
        'Mask and normalized volume dimensions mismatch for subject %s',sid);
    assert(isequal(size(mask),size(vol)),...
        'Mask and volume dimensions mismatch for subject %s',sid);
    assert(isequal(size(mask),size(seg)),...
        'Mask and segmentation volume dimensions mismatch for subject %s',sid);

    %% Apply masks
    try
        fprintf('Applying mask to subject %s\n',sid)
        %%% Apply tissue mask
        % Mask the non-normalized volume (volm)
        if isa(vol,'uint8')
            volm = vol .* uint8(mask);
        elseif isa(vol,'uint16')
            volm = vol .* uint16(mask);
        else
            volm = vol .* mask;
        end
        % Mask the normalized volume (volnm)
        volnm = voln .* uint16(mask);
        % Mask the segmentation (segm)
        segm = logical(seg .* mask);
        
        %%% Clear variables to save memory
        vol = [];
        voln = [];
        seg = [];

        %%% Apply white matter mask to the tissue-masked volumes
        if isa(volm,'uint8')
            volm_wm = volm .* uint8(wm_mask);
        elseif isa(volm,'uint16')
            volm_wm = volm .* uint16(wm_mask);
        else
            volm_wm = volm .* wm_mask;
        end
        volnm_wm = volnm .* uint16(wm_mask);
        segm_wm = logical(segm .* wm_mask);

        %%% Apply gray matter mask (inverse of white matter) to tissue-masked volumes
        if isa(volm,'uint8')
            volm_gm = volm .* uint8(~wm_mask);
        elseif isa(volm,'uint16')
            volm_gm = volm .* uint16(~wm_mask);
        else
            volm_gm = volm .* ~wm_mask;
        end
        volnm_gm = volnm .* uint16(~wm_mask);
        segm_gm = logical(segm .* ~wm_mask);
        fprintf('Finished masking subject %s\n',sid)
    catch
        %%% catch where mask if different data type
        fprintf('failed to apply mask to subject %s\n',sid)
        continue
    end

    
    %% Create overlays
    fprintf('Saving overlays and masked volumes for subject %s\n',sid)
    % Overlay masked normalized volume and segmentation
    fout = fullfile(dpath, subid{ii}, subdir, segdir,...
                    strcat(seg_name, '_refined_mask_overlay_norm.tif'));
    overlay_vol_seg(volnm, segm, 'magenta', fout, false);
    % Overlay masked normalized white matter and segmentation
    fout = fullfile(dpath, subid{ii}, subdir, segdir,...
                    strcat(seg_name, '_wm_refined_mask_overlay_norm.tif'));
    overlay_vol_seg(volnm_wm, segm_wm, 'magenta', fout, false);
    % Overlay masked normalized gray matter and segmentation
    fout = fullfile(dpath, subid{ii}, subdir, segdir,...
                    strcat(seg_name, '_gm_refined_mask_overlay_norm.tif'));
    overlay_vol_seg(volnm_gm, segm_gm, 'magenta', fout, false);
    
    %% Save masked volumes
    %%% Masked non-normalized volume (volm)
    % Entire volume
    volm_out = fullfile(dpath, subid{ii}, mdir,...
        strcat(vol_name,'_refined_masked.tif'));
    segmat2tif(volm, volm_out);
    % White matter
    volwm_out = fullfile(dpath, subid{ii}, mdir,...
        strcat(vol_name,'_wm_refined_masked.tif'));
    segmat2tif(volm_wm, volwm_out);
    % Gray matter
    volgm_out = fullfile(dpath, subid{ii}, mdir,...
        strcat(vol_name,'_gm_refined_masked.tif'));
    segmat2tif(volm_gm, volgm_out);

    %%% Masked normalized volume (volnm)
    vol_out = fullfile(dpath, subid{ii}, mdir,...
        strcat(voln_name,'_refined_masked.tif'));
    segmat2tif(volnm, vol_out);
    vol_out = fullfile(dpath, subid{ii}, mdir,...
        strcat(voln_name,'_wm_refined_masked.tif'));
    segmat2tif(volnm_wm, vol_out);
    vol_out = fullfile(dpath, subid{ii}, mdir,...
        strcat(voln_name,'_gm_refined_masked.tif'));
    segmat2tif(volnm_gm, vol_out);
    
    %% Save masked segmentation
    
    %%% Save the .MAT files
    % Entire volume
    segm_out = fullfile(dpath, subid{ii}, subdir, segdir,...
                        strcat(seg_name,'_refined_masked.mat'));
    save_seg(segm, segm_out);
    % White matter
    segm_out = fullfile(dpath, subid{ii}, subdir, segdir,...
                        strcat(seg_name,'_wm_refined_masked.mat'));
    save_seg(segm_wm, segm_out);
    % Gray matter
    segm_out = fullfile(dpath, subid{ii}, subdir, segdir,...
                        strcat(seg_name,'_gm_refined_masked.mat'));
    save_seg(segm_gm, segm_out);

    %%% Convert to uint8 and save as NIFTI
    % Entire volume
    segm_out = fullfile(dpath, subid{ii}, subdir, segdir,...
                        strcat(seg_name,'_refined_masked.tif'));
    segm_uint8 = uint8(rescale(segm,0,255));
    segmat2tif(segm_uint8, segm_out);
    % White matter
    segm_out = fullfile(dpath, subid{ii}, subdir, segdir,...
                        strcat(seg_name,'_wm_refined_masked.tif'));
    segm_wm_uint8 = uint8(rescale(segm_wm,0,255));
    segmat2tif(segm_wm_uint8, segm_out);
    % Gray matter
    segm_out = fullfile(dpath, subid{ii}, subdir, segdir,...
                        strcat(seg_name,'_gm_refined_masked.tif'));
    segm_gm_uint8 = uint8(rescale(segm_gm,0,255));
    segmat2tif(segm_gm_uint8, segm_out);
    fprintf('\n---------Finished Subject %s---------\n',subid{ii})
end

fprintf('---------Finished Masking all Subjects---------\n')

%% Convert segmentation to graph
subid = {'CTE_7019','CTE_7126','CTE_8572',...
         'NC_8653','NC_21499', 'NC_301181','NC_6974'};
for ii = 1:length(subid)
    fprintf('---------Graphing Subject %s---------\n',subid{ii})
    fprintf('Generating graph data for subject %s\n',subid{ii})
    % Voxel dimensions
    vox_dim = [12, 12, 15];
    % Full filepath to save output
    fullpath = fullfile(dpath, subid{ii}, subdir, segdir);
    % Boolean to view the loop removal debugging plots
    viz = false;
    % Boolean to remove the loops
    rmloop_bool = true;
    
    %%% Entire volume
    % Filename
    fname_seg = strcat(seg_name,'_refined_masked');
    % Import segmentation
    segm = load(fullfile(fullpath, strcat(fname_seg,'.mat')));
    segm = segm.seg;
    % Function to initialize graph and remove loops
    seg_graph_init(segm, vox_dim, fullpath, fname_seg, viz, rmloop_bool)

    %%% White Matter
    % Filename
    fname_seg = strcat(seg_name,'_wm_refined_masked');
    % Import segmentation
    segm_wm = load(fullfile(fullpath, strcat(fname_seg,'.mat')));
    segm_wm = segm_wm.seg;
    % Function to initialize graph and remove loops
    seg_graph_init(segm_wm, vox_dim, fullpath, fname_seg, viz, rmloop_bool)

    %%% Gray Matter
    % Filename
    % Import segmentation
    segm_gm = load(fullfile(fullpath, strcat(fname_seg,'.mat')));
    segm_gm = segm_gm.seg;
    fname_seg = strcat(seg_name,'_gm_refined_masked');
    % Function to initialize graph and remove loops
    seg_graph_init(segm_gm, vox_dim, fullpath, fname_seg, viz, rmloop_bool)

    %%% Debugging Info
    fprintf('---------Finished Subject %s---------\n',subid{ii})
end

%% Function to save segmentation during parallelization
function save_seg(seg, fout)
% save_seg Save the segmentation
% INPUTS:
%   fout (string): output filepath and filename
%   seg (string):  segmentation matrix

save(fout,'seg','-v7.3')

end