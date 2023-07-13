function [ov] = overlay_vol_seg(vol, seg, color)
%overlay_vol_seg Overlay the volume with the segmentation output
% INPUTS:
%   vol (mat): tissue volume (grayscale)
%   seg (mat): segmentations of vasculature (binary)
%   color (string): color of segmentation in overlaid image
% OUTPUTS:
%   ov (mat): segmentations overlaid on top of the tissue volume. The
%           tissue volume will appear grey and the segmentation will appear
%           the color as the input.

%% Overlay each slice in z stack
% Find number of slices in stack
zdim = size(vol, 3);

% Preallocate space for overlaid image
ov = zeros([size(vol), 3]);

% Overlay each slice separately
for ii = 1:zdim
    ov(:,:,ii,:) = imoverlay(vol(:,:,ii), seg(:,:,ii), color);
end

