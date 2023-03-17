function [masked] = apply_mask(orig, mask, fname)
%apply_mask apply 3D mask to image with background
% INPUTS:
%   orig (double matrix): 3D matrix of image with background
%   mask (double matrix): 3D matrix of mask to remove background
%   fname (string): '[output path]\[filename]'
%
% This code was adapted from the Mathworks forum:
% https://www.mathworks.com/matlabcentral/answers/630163-how-can-i-create-
% a-3d-tiff-image-2d-stack-from-a-3d-matrix

%% TODO:
% - run with Hui_Frangi_dataset\200726PSOCT
% - convert TIF files to .MAT
% - write logic to apply mask
% - test on SCC

x = length(orig, 1);
y = length(orig, 2);
z = length(orig, 3);
masked = zeros(x,y,z);

for n = 1:size(orig, 3)    
    % Apply mask
    
    % Generate tiff stack:
    if n == 1
        % First slice:
        imwrite(seg(:, :, n),fname)
    else
        % Subsequent slices:
        imwrite(seg(:, :, n),fname,'WriteMode','append');
    end 
    disp(n)
end

end