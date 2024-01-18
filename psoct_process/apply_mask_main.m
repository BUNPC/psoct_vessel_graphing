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

%% Lookup table of pia pixel intensity reference values
% subid = {'AD_10382', 'AD_20832', 'AD_20969',...
%          'AD_21354', 'AD_21424',...
%          'CTE_6489', 'CTE_6912',...
%          'CTE_7019', 'CTE_8572','CTE_7126',...
%          'NC_6839',  'NC_6974', 'NC_8653',...
%          'NC_21499', 'NC_301181'};

% Matrices that require truncating
% subid = {'CTE_8572','CTE_7126',...
%          'NC_6839',  'NC_6974', 'NC_8653',...
%          'NC_21499', 'NC_301181'};

% Matrices that don't require truncating
subid = {'AD_10382', 'AD_20832', 'AD_20969',...
         'AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912',...
         'CTE_7019'};

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
segdir = 'gsigma_3-5-7_5-7-9_7-9-11';

%%% Filenames
% Combined segmentations (.MAT)
seg_name = 'seg';
% non-normalized PSOCT volume
vol_name = 'ref';
% Normalized PSOCT tissue volume
voln_name = 'ref_4ds_norm_inv';

%% Parallel Processing
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

%% Iterate through subjects. Create and apply mask.
parfor (ii = 1:length(subid), NSLOTS)
% for ii = 1:length(subid)
    %%% Debugging information
    fprintf('Starting Subject %s\n',subid{ii})

    %%% Import volumes, mask, and segmentation
    % Import normalized volume
    fpath = fullfile(dpath, subid{ii}, subdir, strcat(voln_name,'.tif'));
    voln = TIFF2MAT(fpath);   
    
    % Import non-normalized volume
    fpath = fullfile(dpath, subid{ii}, subdir1, strcat(vol_name,'.mat'));
    vol = load(fpath); 
    vol = vol.vol;
    
    % Import mask matrix
    fpath = fullfile(dpath, subid{ii}, subdir1, 'maskec.mat');
    mask = load(fpath);
    mask = mask.mask;
    % Convert mask to logical for matrix operations
    mask = imbinarize(mask);

    % Import combined segmentation file
    fpath = fullfile(dpath, subid{ii}, subdir, subdir2, segdir, 'seg.mat');
    seg = load(fpath,'seg');
    seg = seg.seg;

    %%% Verify dimensions of mask and volumes
    if any(size(mask) ~= size(vol)) || any(size(mask) ~= size(voln))
        % Print debugging info
        fprintf('Dimension mismatch for subject %s\n',subid{ii});
        fprintf('Mask dimensions: [%i,%i,%i]\n',...
            size(mask,1),size(mask,2),size(mask,3));
        fprintf('Volume dimensions: [%i,%i,%i]\n',...
            size(vol,1),size(vol,2),size(vol,3));
        fprintf('Normalized volume dimensions: [%i,%i,%i]\n',...
            size(voln,1),size(voln,2),size(voln,3));
        % Skip this iteration
        continue
    end

    %%% Apply mask
    % Mask the non-normalized volume (volm)
    volm = vol .* mask;
    % Mask the normalized volume (volnm)
    volnm = voln .* uint16(mask);
    % Mask the segmentation (segm)
    segm = seg .* uint8(mask);
    
    %%% Overlays
    % Overlay the segmentation and non-normalized masked volume
    fout = fullfile(dpath, subid{ii}, subdir, subdir2, segdir,...
                    strcat(seg_name, '_refined_mask_overlay.tif'));
    overlay_vol_seg(volm, segm, 'green', fout);

    % Overlay the segmentation and normalized masked volume
    fout = fullfile(dpath, subid{ii}, subdir, subdir2, segdir,...
                    strcat(seg_name, '_refined_mask_overlay_norm.tif'));
    overlay_vol_seg(volnm, segm, 'green', fout);
    
    %%% Save masked volumes
    % Export masked non-normalized volume (volm)
    volm_out = fullfile(dpath, subid{ii}, subdir1,...
        strcat(vol_name,'_refined_masked.tif'));
    save_vol(volm_out, volm);
    
    % Export masked normalized volume (volnm)
    volnm_out = fullfile(dpath, subid{ii}, subdir1,...
        strcat(voln_name,'_refined_masked.tif'));
    save_vol(volnm_out, volnm);
    
    % Export masked segmentation (segm)
    segm_out = fullfile(dpath, subid{ii}, subdir, subdir2, segdir,...
        strcat(seg_name,'_refined_masked.mat'));
    save_seg(segm_out, segm);
    
    %%% Convert to graph
    fprintf('Generating graph data for sub %s\n',subid{ii})
    vox_dim = [12, 12, 15];
    fullpath = fullfile(dpath, subid{ii}, subdir, subdir2, segdir);
    fname_seg = strcat(seg_name,'_refined_masked');
    viz = false;
    rmloop_bool = false;
    seg_graph_init(segm, vox_dim, fullpath, fname_seg, viz, rmloop_bool)

    %%% Debugging Info
    fprintf('Finished Subject %s\n',subid{ii})
end



%% Function to save during parallelization
function save_vol(fout, vol)
% Save the stacked ref matrix
% INPUTS:
%   fout (string): the filepath for the output file to save
%   vol (double matrix): the stack of images to save

% Save the output
save(fout,'vol','-v7.3')

end

%% Function to save during parallelization
function save_seg(fout, seg)
% Save the stacked ref matrix
% INPUTS:
%   fout (string): the filepath for the output file to save
%   vol (double matrix): the stack of images to save

% Save the output
save(fout,'seg','-v7.3')

end

