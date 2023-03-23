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
[orig_re, r, c, s] = reshape_mat(orig);
[mask_re, ~, ~, ~] = reshape_mat(mask);

% Binarize the mask ([0, 255] to [0, 1])
mask_re = uint16(imbinarize(mask_re, 0));

% Element-wise multiplication of mask and volume
masked = uint16(orig_re) .* mask_re;

% Reshape the masked image back to original shape
masked = reshape(masked, [r, c, s]);


%% Reshape matrices (slices appended to rows)
    function [imre, r, c, s] = reshape_mat(im)
        % Concatenate slices to columns.
        % INPUTS:
        %   im (double matrix): 3D matrix to be reshaped
        % OUTPUTS:
        %   imre (double matrix): 2D matrix of reshaped 3D volume
        %   r (double): # rows
        %   c (double): # columns
        %   s (double): # slices
        %
        % Find number rows, columns of volume.
        r = size(im,1);   % number rows
        c =  size(im,2);  % number columns
        s = size(im,3);   % number of slices per volume
        
        % Calculate number of columns for reshaped matrix. Multiply the number of
        % columns per slice by the total number of slices.
        c2 = c .* s;        % (# columns / slice) .* (# slices)
        
        % Reshape the volume to have the number of columns as the c2
        imre = reshape(im, [r,c2]);
    end

end