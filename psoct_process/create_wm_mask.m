function [wm_mask, gm_mask, tissue_mask] = create_wm_mask(mus,oct,th,viz)
%CREATE_WM_MASK create white matter mask from scattering map
%
% INPUTS:
%   mus (matrix): scattering map
%   oct (matrix): OCT intensity matrix
%   th (matrix): min threshold for segmenting white matter in scatter map
%   viz (bool): visualize the debugging plots
% OUTPUTS:
%   wm_mask (matrix): 
%   gm_mask (matrix): 
%   tissue_mask (matrix): 

%% Elements for modifying mask
se_erode = strel('disk', 10);
tiss_erode = strel('disk', 20);

%% Segment WM from mus map        
% Threshold and take the largest segment
mus_th = mus;
mus_th(mus_th<=th) = 0;
mus_th = bwareaopen(mus_th, 50);

% Fill the holes
close_mus = imfill(mus_th,'holes');

% Erode the border and remove voxel groups smaller than nmin
wm_mask = imerode(close_mus,se_erode);
nmin = 10000;
wm_mask = bwareaopen(wm_mask, nmin);

%%% Debugging figures
if viz
    figure; imshow(mus_th); title('Largest Segment');
    figure; imshow(close_mus); title('Close holes')
    figure; imshow(wm_mask); title('WM Mask')
    % Overlay the mask with the OCT
    figure; h = imshow(oct);
    set(h, 'AlphaData',wm_mask);
end
%% Segment GM from inverse of WM mask
% Create tissue mask from scattering mask, given the knowledge that the
% agarose scattering coeff. is approximately 0.
tissue_mask = mus;
tissue_mask(tissue_mask > 0.01) = 1;

% Erode the borders
tissue_mask = imerode(tissue_mask, tiss_erode);

% Close the holes
tissue_mask = imfill(tissue_mask, 'holes');

% Gray matter mask. Inverse of white matter mask and tissue mask
gm_mask = bsxfun(@times, tissue_mask, cast(~wm_mask,'like',tissue_mask));
gm_mask = logical(gm_mask);

% Remove voxels with fewer than nmin connections
gm_mask = bwareaopen(gm_mask, nmin);

%%% Manually modify the mask
% immask(gm_mask);

%%% Debugging Figures
if viz
    figure; imagesc(tissue_mask); title('Tissue Mask');
    figure; imshow(gm_mask); title('GM mask');
    % Close all figures
    close all
end
end