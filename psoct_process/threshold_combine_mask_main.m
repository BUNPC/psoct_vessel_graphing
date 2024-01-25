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
    %%% Init paralell pool
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
% Subfolders for sigma arrays
sigmas = {'gsigma_1-3-5_gsize_5-13-21','gsigma_2-3-4_gsize_9-13-17',...
        'gsigma_3-5-7_gsize_13-21-29','gsigma_5-7-9_gsize_21-29-37',...
        'gsigma_7--9-11_gsize_29-37-45'};
sigma_field = {'sigma135','sigma234','sigma357','sigma579','sigma7911'};
sigout = 'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11';

%%% Filenames
% Probability matrix (.MAT)
pname = 'probability_map.mat';
% Normalized PSOCT tissue volume
voln_name = 'ref_4ds_norm_inv';

%%% Array of probability thresholds for segmentation
th = [0.18, 0.20, 0.22, 0.24];
pfield = {'p18','p20','p22','p24'};

%% Iterate subjects, threhsold prob, combine segs, apply mask.
parfor (ii = 1:length(subid), NSLOTS)
    %%% Initialize variables for parfor
    segs = struct();
    pmat = [];
    
    %%% Debugging information
    fprintf('---------Starting Subject %s---------\n',subid{ii})
    
    %%% Iterate over sigma arrays
    fprintf('thresholding probability matrices')
    for j = 1:length(sigmas)
        % Load the probability map
        fpath = fullfile(dpath, subid{ii}, subdir, sigmas{j}, pname);
        pmat = load(fpath,'pmat');
        % Convert to single (save memory)
        pmat = single(pmat.pmat);
        
        %%% Iterate over probability thresholds
        for k = 1:length(th)
            % Apply threshold
            pth = pmat;
            pth(pth < th(k)) = 0;
            pth(pth >= th(k)) = 1;
            % Save to struct
            segs.(pfield{k}).(sigma_field{j}) = uint8(pth);
        end
    end

    %%% Iterate probability threshold and combine segments across sigmas
    fprintf('combining segmentations')
    for j = 1:length(th)
        % Init matrix to store segmentations
        segmat = uint8(zeros(size(pmat)));
        % Iterate across sigma array segmentations        
        for k = 1:length(sigmas)
            segmat = segmat +...
                segs.(pfield{j}).(sigma_field{k});
        end
        % Add combined segmentation to struct
        segs.(pfield{j}).combined = segmat;
    end
    
    %% Iterate threhold, apply mask to combined segmentations
    %%% Import volumes, mask, and segmentation
    % Import normalized volume
    fprintf('importing normalized volume\n')
    fpath = fullfile(dpath, subid{ii}, subdir, strcat(voln_name,'.tif'));
    voln = TIFF2MAT(fpath);   
    % Import mask matrix
    fprintf('importing mask\n')
    fpath = fullfile(dpath, subid{ii}, subdir1, 'maskec.mat');
    mask = load(fpath);
    mask = mask.mask;
    % Convert mask to logical for matrix operations
    mask = imbinarize(mask);

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
    % The segmentation file was made using the normalized volume. The mask
    % was created from the ref[#].mat files, and some have one extra voxel
    % than the normalied file. Therefore, the dimensions of the mask must
    % be made equal to those of the volume.
    if any(size(mask) ~= size(voln))
        %%% Truncate mask & vol (ref) to match dimensions of normalized
        m = size(mask);
        vn = size(voln);
        % Find smallest dimensions
        ymin = min(m(1), vn(1));
        xmin = min(m(2), vn(2));
        zmin = min(m(3), vn(3));
        % Truncate all volumes
        mask = mask(1:ymin, 1:xmin, 1:zmin);
        voln = voln(1:ymin, 1:xmin, 1:zmin);
        % Assert that the new matrices are all equivalent
        assert(isequal(size(mask), size(voln)),...
            'Mask and voln dimensions mismatch')
    end
    
    %%% Apply mask to normalized volume
    fprintf('apply mask\n')
    if isa(voln, 'uint16')
        % Mask the normalized volume (volnm)
        volnm = voln .* uint16(mask);   
    elseif isa(voln, 'uint8')
        % Mask the normalized volume (volnm)
        volnm = voln .* uint8(mask);
    else
        fprintf('unexpected datatype of voln for sub %s',subid{ii})
    end
    %%% Export masked normalized volume (volnm)
    fprintf('exporting masked normalized volume\n')
    fout = fullfile(dpath,subid{ii},subdir,...
            strcat(voln_name,'_refined_masked.tif'));
    segmat2tif(volnm, fout);
    
    %%% Iterate over thresholds
    for j = 1:length(th)
        % Combined segmentation matrix
        seg = segs.(pfield{j}).combined;
        % Truncate dimensions of segmentation to match mask and voln
        if any(size(seg) ~= size(mask))        
            seg = seg(1:ymin, 1:xmin, 1:zmin);
            assert(isequal(size(seg), size(voln)),...
            'Seg and voln dimensions mismatch')
        end      
    
        %% Apply mask to segmentation
        % Mask the segmentation (segm)
        segm = seg .* uint8(mask);
        
        %% Create overlays, save outputs, create/save graph
        %%% Overlays the segmentation and normalized masked volume
        fprintf('overlaying segmentation and masked normalized volume\n')
        fout = fullfile(dpath,subid{ii},subdir,'combined_segs',sigout,...
            pfield{j}, 'seg_refined_mask_overlay_norm.tif');
        overlay_vol_seg(volnm, segm, 'green', fout);
              
        %%% Export masked segmentation (segm)
        fprintf('exporting masked segmentation\n')
        fout = fullfile(dpath,subid{ii},subdir,'combined_segs',sigout,...
            pfield{j},strcat('seg_refined_masked.mat'));
        save_vol(fout,segm);
        fout = fullfile(dpath,subid{ii},subdir,'combined_segs',sigout,...
            pfield{j},strcat('seg_refined_masked.tif'));
        segmat2tif(segm, fout);
        
        %%% Convert to graph
        fprintf('Generating graph data for sub %s\n',subid{ii})
        vox_dim = [12, 12, 15];
        fout = fullfile(dpath,subid{ii},subdir,'combined_segs',sigout,...
            pfield{j});
        fname_seg = 'seg_refined_masked';
        viz = false;
        rmloop_bool = false;
        seg_graph_init(segm, vox_dim, fout, fname_seg, viz, rmloop_bool)
    end
    
    %%% Debugging Info
    fprintf('---------Finished Subject %s---------\n',subid{ii})
end

function save_vol(fout, vol)
% Save volume or segmentation
% INPUTS:
%   fout (string): filepath output
%   vol (matrix): segmentation or volume to save
save(fout, 'vol', '-v7.3');

end
