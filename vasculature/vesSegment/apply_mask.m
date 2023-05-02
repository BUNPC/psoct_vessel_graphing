function [masked] = apply_mask(orig, mask)
%apply_mask apply 3D mask to image with background
% INPUTS:
%   orig (double matrix): 3D matrix of image with background
%   mask (double matrix): 3D matrix of mask to remove background
% OUTPUTS:
%   masked (double matrix): 2D matrix of masked original image
%
% This code was adapted from the Mathworks forum:
% https://www.mathworks.com/matlabcentral/answers/630163-how-can-i-create-
% a-3d-tiff-image-2d-stack-from-a-3d-matrix

%% Reshape 3D --> 2D, multiply, reshape 2D --> 3D

% Reshape mask and volume
[orig_re, r, c, s] = reshape_3dmat_to_2d(orig);
[mask_re, ~, ~, ~] = reshape_3dmat_to_2d(mask);

% Binarize the mask ([0, 255] to [0, 1])
mask_re = uint16(imbinarize(mask_re, 0));

% Element-wise multiplication of mask and volume
masked = uint16(orig_re) .* mask_re;

% Reshape the masked image back to original shape
masked = reshape(masked, [r, c, s]);

end