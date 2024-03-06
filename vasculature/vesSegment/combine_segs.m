function combine_segs(ov)
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
%%% OUTPUTS:
%       segfinal (matrix): combined segmentation
%       graphfinal (struct): graph from combined segmentation

% TODO: 
% - pass in the volume .TIF file
% - create output file path

%% Overlay all segmentations
% Iterate over subjects
for ii = 1:size(ov,2)
    % Load segmentations from first sigma array
    seg = TIFF2MAT(ov(ii).f(1).fpath);

    % Iterate over segmentations (start at second)
    for j = 2:length(ov(ii).s)
        % Load segmentation into local variable
        tmp = TIFF2MAT(ov(ii).f(j).fpath);
        % Overlay the segmentations
        seg = seg + tmp;
    end
    
    %%% Save segmentation output
    % Create folder if it doesn't exist
    if ~exist(ov(ii).dirout,'dir')
        mkdir(ov(ii).dirout);
    end
    % Filename
    segout = 'seg.mat';
    segout = fullfile(ov(ii).dirout, segout);
    % Save combined segmentations
    save(segout,'seg','-v7.3')


    %%% Generate quality assurance overlay with volume
    % Import volume
    vol = TIFF2MAT(ov(ii).vol);
    % Select color for segmentation
    color = 'green';
    % Output file path
    fout = fullfile(ov(ii).dirout, 'segs_combined_overlay.tif');
    % Call function to overlay segmentation with volume and save output
    overlay_vol_seg(vol, seg, color, fout)

end

end