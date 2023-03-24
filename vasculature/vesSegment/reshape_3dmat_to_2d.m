function [imre, r, c, s] = reshape_3dmat_to_2d(im)
%RESHAPE_3DMAT_TO_2D Concatenate slices to columns.
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

