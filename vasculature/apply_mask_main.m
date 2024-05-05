%% Apply the tissue, white matter, and gray matter masks to segmentation
% This script applies several masks to the segmentation. These are the
% following masks: tissue, white matter (WM), gray matter (GM), sulci, and
% gyri. The names of the masks are: mask_tiss, mask_gm, mask_wm,
% mask_sulci, mask_gyri. The purpose of applying these masks is to analyze
% the various regional characteristics of the brain tissue specimens.
% The masks are all stored in the directory:
% Ann_Mckee_samples_55T/[subID]/dist_corrected/volume/ref/masks/
%
% Then, this script converts the masked components of segmentation into
% graphs before analyzing the graphs.

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
topdir = mydir(1:idcs(end));
addpath(genpath(topdir));

%% Initialize paralell pool
% Set # threads = # cores for job
NSLOTS = str2num(getenv('NSLOTS'));
maxNumCompThreads(NSLOTS);
% Check to see if we already have a parpool, if not create one with
% our desired parameters
poolobj = gcp('nocreate');
if isempty(poolobj)
    myCluster=parcluster('local');
    % Ensure multiple parpool jobs don't overwrite other's temp files
    myCluster.JobStorageLocation = getenv('TMPDIR');
    poolobj = parpool(myCluster, NSLOTS);
end

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

%%% Minimum number of voxels to retain a segmentations
vox_min = 100;

%% Initialize directories and filenames

%%% Directories 
% Upper level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T';    
% Subfolder with normalized volume
subdir = '/dist_corrected/volume';
% Subfolder containing non-normalized ref files
subdir1 = '/dist_corrected/volume/ref';
% Combined segmentation subfolder
segdir = '/combined_segs/gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/';
segdir2 = append('vox_min_',num2str(vox_min));
% Mask subfolder
mdir = '/dist_corrected/volume/ref/masks';

%%% Filenames
% Combined segmentations (.MAT)
seg_name = 'seg';
% non-normalized PSOCT volume
vol_name = 'ref';
% Normalized PSOCT tissue volume
voln_name = 'ref_4ds_norm_inv';

%% Iterate through subjects. Create and apply mask_tiss.
parfor (ii = 3:length(subid),NSLOTS)
% for ii = 2:length(subid)
    %%% Debugging information
    fprintf('\n---------Starting Subject %s---------\n',subid{ii})

    %% Import tissue volumes, segmentation, and masks
    %%% Import normalized volume
    fpath = fullfile(dpath, subid{ii}, subdir, strcat(voln_name,'.tif'));
    voln = TIFF2MAT(fpath);   

    %%% Import non-normalized volume
    fpath = fullfile(dpath, subid{ii}, subdir1, strcat(vol_name,'.mat'));
    vol = load(fpath); 
    vol = vol.vol;

    %%% Import combined segmentation file
    fpath = fullfile(dpath, subid{ii}, subdir, segdir, 'seg.mat');
    seg = import_volume(fpath);
    % Remove the components with fewer than 100 voxels in connectivity
    seg = rm_short_vessels(seg, vox_min);

    %%% Import tissue mask
    fpath = fullfile(dpath, subid{ii}, mdir, 'mask_tiss.mat');
    mask_tiss = import_volume(fpath);

    %%% Import white matter mask
    fpath = fullfile(dpath, subid{ii}, mdir,'mask_wm.mat');
    mask_wm = import_volume(fpath);
    % Bitwise operation between WM and tissue masks.
    mask_wm = bsxfun(@times, mask_tiss, cast(mask_wm,'like',mask_tiss));
    mask_wm = logical(mask_wm);

    %%% Import gray matter mask
    fpath = fullfile(dpath, subid{ii}, mdir,'mask_gm.mat');
    mask_gm = import_volume(fpath);
    % Bitwise operation between GM and tissue masks.
    mask_gm = bsxfun(@times, mask_tiss, cast(mask_gm,'like',mask_tiss));
    mask_gm = logical(mask_gm);

    %%% Import sulci & gyri masks
    fpath = fullfile(dpath, subid{ii}, mdir,'mask_sulci.mat');
    mask_sulci = import_volume(fpath);
    fpath = fullfile(dpath, subid{ii}, mdir,'mask_gyri.mat');
    mask_gyri = import_volume(fpath);
    
    %% Verify dimensions of masks and volumes
    % The segmentation file was made using the normalized volume. Therefore
    % the dimensions of the mask_tiss must be made equal to those of the seg.

    %%% Set zmin for truncating finals slices
    % Truncate the mask_tiss, and the others will be subsequently truncated.
    sid = subid{ii};
    % Retrieve the zmin from the array
    zmin = d({sid});
    % Truncate the tissue mask_tiss
    mask_tiss = mask_tiss(:,:,1:zmin);

    %%% Truncate mask_tiss & vol (ref) to match dimensions of normalized
    % Find smallest dimensions
    ymin = min([size(mask_tiss,1), size(mask_wm,1), size(voln,1)]);
    xmin = min([size(mask_tiss,2), size(mask_wm,2), size(voln,2)]);
    zmin = min([size(mask_tiss,3), size(mask_wm,3), size(voln,3)]);
    % Truncate all volumes
    mask_tiss =     mask_tiss(1:ymin, 1:xmin, 1:zmin);
    mask_wm =       mask_wm(1:ymin, 1:xmin, 1:zmin);
    mask_gm =       mask_gm(1:ymin, 1:xmin, 1:zmin);
    mask_sulci =    mask_sulci(1:ymin, 1:xmin, 1:zmin);
    mask_gyri =     mask_gyri(1:ymin, 1:xmin, 1:zmin);
    vol =           vol( 1:ymin, 1:xmin, 1:zmin);
    voln =          voln(1:ymin, 1:xmin, 1:zmin);
    seg =           seg( 1:ymin, 1:xmin, 1:zmin);
    % Assert that the new matrices are all equivalent
    assert(isequal(size(mask_tiss),size(voln)),...
        'Mask and normalized volume dimensions mismatch for subject %s',sid);
    assert(isequal(size(mask_tiss),size(vol)),...
        'Mask and volume dimensions mismatch for subject %s',sid);
    assert(isequal(size(mask_tiss),size(seg)),...
        'Mask and segmentation volume dimensions mismatch for subject %s',sid);

    %% Apply masks and save outputs
    try
        %% Apply masks
        % Meaning of variables:
        %   vol = OCT volume
        %   voln = normalized OCT volume
        %   volnm = normalized OCT volume (tissue masked)
        %   volnm_wm = normalized OCT volume (tissue, WM masked)
        %   volnm_gm = normalized OCT volume (tissue, GM masked)
        %   seg = segmentation (not masked)
        %   segm = segmentation (tissue masked)
        %   segm_wm = segmentation (tissue, WM masked)
        %   segm_gm = segmentation (tissue, GM masked)

        fprintf('Applying masks to subject %s\n',sid)
        %%% Apply tissue mask
        % Mask the non-normalized volume (volm)
        if isa(vol,'uint8')
            volm = vol .* uint8(mask_tiss);
        elseif isa(vol,'uint16')
            volm = vol .* uint16(mask_tiss);
        else
            volm = vol .* mask_tiss;
        end
        % Mask the normalized volume (volnm)
        volnm = voln .* uint16(mask_tiss);
        % Mask the segmentation (segm)
        segm = logical(seg .* mask_tiss);
        
        %%% Clear variables to save memory
        vol = [];
        voln = [];
        seg = [];

        %%% Apply masks (WM, GM, sulci, gyri) to non-normalized OCT
        if isa(volm,'uint8')
            volm_wm = volm .* uint8(mask_wm);
            volm_gm = volm .* uint8(mask_gm);
            volm_sulci = volm .* uint8(mask_sulci);
            volm_gyri = volm .* uint8(mask_gyri);
        elseif isa(volm,'uint16')
            volm_wm = volm .* uint16(mask_wm);
            volm_gm = volm .* uint16(mask_gm);
            volm_sulci = volm .* uint16(mask_sulci);
            volm_gyri = volm .* uint16(mask_gyri);
        else
            volm_wm = volm .* mask_wm;
            volm_gm = volm .* mask_gm;
            volm_sulci = volm .* mask_sulci;
            volm_gyri = volm .* mask_gyri;
        end

        %%% Apply masks (WM,GM,sulci,gyri) to normalized OCT
        volnm_wm =          volnm .*        uint16(mask_wm);
        volnm_gm =          volnm .*        uint16(mask_gm);
        volnm_sulci =       volnm .*        uint16(mask_sulci);
        volnm_gyri =        volnm .*        uint16(mask_gyri);
        volnm_wm_sulci =    volnm_sulci .*  uint16(mask_wm);
        volnm_wm_gyri =     volnm_gyri .*   uint16(mask_wm);
        volnm_gm_sulci =    volnm_sulci .*  uint16(mask_gm);
        volnm_gm_gyri =     volnm_gyri .*   uint16(mask_gm);
        
        %%% Apply masks (WM,GM,sulci,gyri) to segmentation
        segm_wm =       logical(segm .*         mask_wm);
        segm_gm =       logical(segm .*         mask_gm);
        segm_sulci =    logical(segm .*         mask_sulci);
        segm_gyri =     logical(segm .*         mask_gyri);
        segm_wm_sulci = logical(segm_sulci .*   mask_wm);
        segm_wm_gyri =  logical(segm_gyri .*    mask_wm);
        segm_gm_sulci = logical(segm_sulci .*   mask_gm);
        segm_gm_gyri =  logical(segm_gyri .*    mask_gm);

        fprintf('Finished masking subject %s\n',sid)
        %% Create overlays
        fprintf('Saving overlays and masked volumes for subject %s\n',sid)

        % Create the full path to output overlays
        masked_seg_output = fullfile(dpath,subid{ii},subdir,segdir,segdir2);
        if ~exist(masked_seg_output,'dir')
            mkdir(masked_seg_output)
        end

        % Overlay masked normalized volume and segmentation
        fout = fullfile(masked_seg_output,...
                strcat(seg_name, '_refined_mask_overlay_norm.tif'));
        overlay_vol_seg(volnm, segm, 'magenta', fout, false);
        % Overlay masked normalized WM and segmentation
        fout = fullfile(masked_seg_output,...
                strcat(seg_name, '_wm_refined_mask_overlay_norm.tif'));
        overlay_vol_seg(volnm_wm, segm_wm, 'magenta', fout, false);
        % Overlay masked normalized GM and segmentation
        fout = fullfile(masked_seg_output,...
                strcat(seg_name, '_gm_refined_mask_overlay_norm.tif'));
        overlay_vol_seg(volnm_gm, segm_gm, 'magenta', fout, false);
        % Overlay SULCI and segmentation
        fout = fullfile(masked_seg_output,...
                strcat(seg_name, '_sulci_refined_mask_overlay_norm.tif'));
        overlay_vol_seg(volnm_sulci, segm_sulci, 'magenta', fout, false);
        % Overlay GYRI and segmentation
        fout = fullfile(masked_seg_output,...
                strcat(seg_name, '_gyri_refined_mask_overlay_norm.tif'));
        overlay_vol_seg(volnm_gyri, segm_gyri, 'magenta', fout, false);
        
        %% Save masked volumes (without segmentation)
        %{
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
        % Sulci
        vol_sulci_out = fullfile(dpath, subid{ii}, mdir,...
            strcat(vol_name,'_sulci_refined_masked.tif'));
        segmat2tif(volm_sulci, vol_sulci_out);
        % Gyri
        vol_gyri_out = fullfile(dpath, subid{ii}, mdir,...
            strcat(vol_name,'_gyri_refined_masked.tif'));
        segmat2tif(volm_gyri, vol_gyri_out);
    
        %%% Masked normalized volume (volnm)
        % Tissue Mask
        vol_out = fullfile(dpath, subid{ii}, mdir,...
            strcat(voln_name,'_refined_masked.tif'));
        segmat2tif(volnm, vol_out);
        % White Matter
        vol_out = fullfile(dpath, subid{ii}, mdir,...
            strcat(voln_name,'_wm_refined_masked.tif'));
        segmat2tif(volnm_wm, vol_out);
        % Gray Matter
        vol_out = fullfile(dpath, subid{ii}, mdir,...
            strcat(voln_name,'_gm_refined_masked.tif'));
        segmat2tif(volnm_gm, vol_out);
        % Sulci
        vol_sulci_out = fullfile(dpath, subid{ii}, mdir,...
            strcat(voln_name,'_sulci_refined_masked.tif'));
        segmat2tif(volnm_sulci, vol_sulci_out);
        % Gyri
        vol_gyri_out = fullfile(dpath, subid{ii}, mdir,...
            strcat(voln_name,'_gyri_refined_masked.tif'));
        segmat2tif(volnm_gyri, vol_gyri_out);
        %}
        %% Save masked segmentation as .MAT
        % Create the full path to output overlays
        masked_seg_output = fullfile(dpath,subid{ii},subdir,segdir,segdir2);
        if ~exist(masked_seg_output,'dir')
            mkdir(masked_seg_output)
        end
        
        % Entire volume
        fout = fullfile(masked_seg_output,...
                        strcat(seg_name,'_refined_masked.mat'));
        save_seg(segm, fout);

        % White matter
        fout = fullfile(masked_seg_output,...
                        strcat(seg_name,'_wm_refined_masked.mat'));
        save_seg(segm_wm, fout);

        % Gray matter
        fout = fullfile(masked_seg_output,...
                        strcat(seg_name,'_gm_refined_masked.mat'));
        save_seg(segm_gm, fout);

        % Sulci
        fout = fullfile(masked_seg_output,...
                        strcat(seg_name,'_sulci_refined_masked.mat'));
        save_seg(segm_sulci, fout);

        % Sulci - WM
        fout = fullfile(masked_seg_output,...
                        strcat(seg_name,'_wm_sulci_refined_masked.mat'));
        save_seg(segm_wm_sulci, fout);

        % Sulci - GM
        fout = fullfile(masked_seg_output,...
                        strcat(seg_name,'_gm_sulci_refined_masked.mat'));
        save_seg(segm_gm_sulci, fout);

        % Gyri
        fout = fullfile(masked_seg_output,...
                        strcat(seg_name,'_gyri_refined_masked.mat'));
        save_seg(segm_gyri, fout);

        % Gyri - WM
        fout = fullfile(masked_seg_output,...
                        strcat(seg_name,'_wm_gyri_refined_masked.mat'));
        save_seg(segm_wm_gyri, fout);

        % Gyri - GM
        fout = fullfile(masked_seg_output,...
                        strcat(seg_name,'_gm_gyri_refined_masked.mat'));
        save_seg(segm_gm_gyri, fout);

        %% Convert masked segmentation as TIF
        % Entire volume
        fout = fullfile(masked_seg_output,...
                        strcat(seg_name,'_refined_masked.tif'));
        segmat2tif(uint8(rescale(segm,0,255)), fout);

        % White matter
        fout = fullfile(masked_seg_output,...
                        strcat(seg_name,'_wm_refined_masked.tif'));
        segmat2tif(uint8(rescale(segm_wm,0,255)), fout);

        % Gray matter
        fout = fullfile(masked_seg_output,...
                        strcat(seg_name,'_gm_refined_masked.tif'));
        segmat2tif(uint8(rescale(segm_gm,0,255)), fout);

        % Sulci
        fout = fullfile(masked_seg_output,...
                        strcat(seg_name,'_sulci_refined_masked.tif'));
        segmat2tif(uint8(rescale(segm_sulci,0,255)), fout);

        % Sulci - WM
        fout = fullfile(masked_seg_output,...
                        strcat(seg_name,'_wm_sulci_refined_masked.tif'));
        segmat2tif(uint8(rescale(segm_wm_sulci,0,255)), fout);

        % Sulci - GM
        fout = fullfile(masked_seg_output,...
                        strcat(seg_name,'_gm_sulci_refined_masked.tif'));
        segmat2tif(uint8(rescale(segm_gm_sulci,0,255)), fout);

        % Gyri
        fout = fullfile(masked_seg_output,...
                        strcat(seg_name,'_gyri_refined_masked.tif'));
        segmat2tif(uint8(rescale(segm_gyri,0,255)), fout);

        % Gyri - WM
        fout = fullfile(masked_seg_output,...
                        strcat(seg_name,'_wm_gyri_refined_masked.tif'));
        segmat2tif(uint8(rescale(segm_wm_gyri,0,255)), fout);

        % Gyri - GM
        fout = fullfile(masked_seg_output,...
                        strcat(seg_name,'_gm_gyri_refined_masked.tif'));
        segmat2tif(uint8(rescale(segm_gm_gyri,0,255)), fout);

        fprintf('\n---------Finished Subject %s---------\n',subid{ii})

        %% Clear Variables to save memory
        volnm_wm =          [];
        volnm_gm =          [];
        volnm_sulci =       [];
        volnm_gyri =        [];
        volnm_wm_sulci =    [];
        volnm_wm_gyri =     [];
        volnm_gm_sulci =    [];
        volnm_gm_gyri =     [];
        segm_wm =       [];
        segm_gm =       [];
        segm_sulci =    [];
        segm_gyri =     [];
        segm_wm_sulci = [];
        segm_wm_gyri =  [];
        segm_gm_sulci = [];
        segm_gm_gyri =  [];
    
    catch
        %%% catch where mask failed, if different data type
        fprintf('failed to apply mask to subject %s\n',sid)
        continue
    end
end

fprintf('---------Finished Masking all Subjects---------\n')

%% Convert segmentation to graph
for ii = 1:length(subid)
    fprintf('---------Graphing Subject %s---------\n',subid{ii})
    % Voxel dimensions
    vox_dim = [12, 12, 15];
    % Full filepath to save output
    fullpath = fullfile(dpath, subid{ii}, subdir, segdir, segdir2);
    % Boolean to view the loop removal debugging plots
    viz = false;
    % Boolean to remove the loops
    rmloop_bool = true;
    
    %%% Entire volume
    % Import segmentation
    fname_seg = strcat(seg_name,'_refined_masked');
    seg = import_volume(fullfile(fullpath, strcat(fname_seg,'.mat')));
    % Initialize graph and remove loops
    seg_graph_init(seg, vox_dim, fullpath, fname_seg, viz, rmloop_bool)

    %%% White Matter
    % Import segmentation
    fname_seg = strcat(seg_name,'_wm_refined_masked');
    seg = import_volume(fullfile(fullpath, strcat(fname_seg,'.mat')));
    % Initialize graph and remove loops
    seg_graph_init(seg, vox_dim, fullpath, fname_seg, viz, rmloop_bool)

    %%% Gray Matter
    % Import segmentation
    fname_seg = strcat(seg_name,'_gm_refined_masked');
    seg = import_volume(fullfile(fullpath, strcat(fname_seg,'.mat')));
    % Initialize graph and remove loops
    seg_graph_init(seg, vox_dim, fullpath, fname_seg, viz, rmloop_bool)

    %%% Sulci
    % Import segmentation
    fname_seg = strcat(seg_name,'_sulci_refined_masked');
    seg = import_volume(fullfile(fullpath, strcat(fname_seg,'.mat')));
    % Initialize graph and remove loops
    seg_graph_init(seg, vox_dim, fullpath, fname_seg, viz, rmloop_bool)

    %%% Gyri
    % Import segmentation
    fname_seg = strcat(seg_name,'_gyri_refined_masked');
    seg = import_volume(fullfile(fullpath, strcat(fname_seg,'.mat')));
    % Initialize graph and remove loops
    seg_graph_init(seg, vox_dim, fullpath, fname_seg, viz, rmloop_bool)

    %%% Sulci - WM
    % Import segmentation
    fname_seg = strcat(seg_name,'_wm_sulci_refined_masked');
    seg = import_volume(fullfile(fullpath, strcat(fname_seg,'.mat')));
    % Initialize graph and remove loops
    seg_graph_init(seg, vox_dim, fullpath, fname_seg, viz, rmloop_bool)

    %%% Gyri - WM
    % Import segmentation
    fname_seg = strcat(seg_name,'_wm_gyri_refined_masked');
    seg = import_volume(fullfile(fullpath, strcat(fname_seg,'.mat')));
    % Initialize graph and remove loops
    seg_graph_init(seg, vox_dim, fullpath, fname_seg, viz, rmloop_bool)

    %%% Sulci - GM
    % Import segmentation
    fname_seg = strcat(seg_name,'_gm_sulci_refined_masked');
    seg = import_volume(fullfile(fullpath, strcat(fname_seg,'.mat')));
    % Initialize graph and remove loops
    seg_graph_init(seg, vox_dim, fullpath, fname_seg, viz, rmloop_bool)

    %%% Gyri - GM
    % Import segmentation
    fname_seg = strcat(seg_name,'_gm_gyri_refined_masked');
    seg = import_volume(fullfile(fullpath, strcat(fname_seg,'.mat')));
    % Initialize graph and remove loops
    seg_graph_init(seg, vox_dim, fullpath, fname_seg, viz, rmloop_bool)

    %%% Debugging Info
    fprintf('---------Finished Graphing Subject %s---------\n',subid{ii})
end

%% Function to save segmentation during parallelization
function save_seg(seg, fout)
% save_seg Save the segmentation
% INPUTS:
%   fout (string): output filepath and filename
%   seg (string):  segmentation matrix

save(fout,'seg','-v7.3')

end

%% Function to import masks
function mask = import_volume(fpath)
% import_volume Load the .MAT of the respective mask
% INPUTS:
%   fpath (string): path to the mask file
% OUTPUTS:
%   mask (logical matrix): binary matrix of mask (1=keep, 0=discard)

% Import the mask file
mask = load(fpath);
field = fieldnames(mask);
field = field{:};
mask = mask.(field);

% Convert mask to logical for matrix operations
if ~isa(mask,'logical')
    mask = imbinarize(mask);
end

end