function combine_pmats(ov, th, subdir)
% Combine segmentations for the same volume from multiple seg params.
%   The Frangi vessel enhancement algorithm uses a Gaussian filter to
%   enhance vessels of a specific diameter. It is necessary to run the
%   Frangi with multiple Gaussian kernels to segment a range of vessel
%   diameters. This function will combine the segmentations from each of
%   these segmentations, to generate a segmentation that includes small to
%   large vessels.
%%% INPUTS:
%       ov (struct): struct containing filepath
%           ov(#).subid: subject ID
%           ov(#).basepath: base filepath to segmentation files
%           ov(#).fpath: cell array of strings of full file path to the
%                       .TIFs for each segmentation
%       th (double): probability threshold
%%% OUTPUTS:
%       segfinal (matrix): combined segmentation
%       graphfinal (struct): graph from combined segmentation

%% Init paralell pool
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

%% Overlay all segmentations
% Iterate over subjects
parfor (ii = 1:size(ov,2), NSLOTS)
    % Debugging
    fprintf('----Starting subject %s----\n',ov(ii).subid)
    % Load segmentations from first sigma array
    pmat = load(ov(ii).f(1).fpath);
    pmat = pmat.pmat;
    % Threshold the probability matrix
    seg = pmat_threshold(pmat, th);

    % Iterate over segmentations (start at second)
    for j = 2:length(ov(ii).s)
        % Load next sigma value
        pmat = load(ov(ii).f(j).fpath);
        pmat = pmat.pmat;
        % Threshold and convert to logical
        tmp = pmat_threshold(pmat, th);
        % Overlay the segmentations
        seg = logical(seg + tmp);
    end
    % Convert combined segmentation to logical
    seg = logical(seg);
    
    %%% Save segmentation output
    % Output directory
    dirout = fullfile(ov(ii).dirout,subdir);
    % Create folder if it doesn't exist
    if ~exist(dirout,'dir')
        mkdir(dirout);
    end
    % Filename
    segout = fullfile(dirout, 'seg.mat');
    % Save combined segmentations
    save_seg(segout, seg);

    %%% Generate quality assurance overlay with volume
    % Import volume
    vol = TIFF2MAT(ov(ii).vol);
    % Select color for segmentation
    color = 'magenta';
    % Output file path
    fout = fullfile(dirout, 'segs_combined_overlay.tif');
    % Call function to overlay segmentation with volume and save output
    overlay_vol_seg(vol, seg, color, fout, false)

    fprintf('----Finished subject %s----\n',ov(ii).subid)
end



end
%% Function to threshold probability matrix during parallelization
function [pmat] = pmat_threshold(pmat, th)
    % Threshold the probability matrix
    % INPUTS:
    %       pmat (double matrix): probability matrix
    %       th (double): minimum probability threshold
    % OUTPUTS:
    %       pmat_th (logical matrix): boolean matrix of segmentation
    
    % Threshold the probability matrix
    pmat(pmat<th) = 0;
    pmat(pmat>=th) = 1;
    % Convert probability matrix to logical to save space
    pmat = logical(pmat);
end


%% Function to save during parallelization
function save_seg(fout, seg)
% Save the segmentation
% INPUTS:
%   fout (string): the filepath for the output file
%   seg (double matrix): the segmentation matrix to save

save(fout,'seg','-v7.3')
end
