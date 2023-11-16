function rgb_stack_2_tif(stack, fname)
%rgb_stack_2_tif convert 3D matrix of segmentation to a 2D stack TIF file
% INPUTS:
%   stack (double matrix): 4D stack of RGB images
%   fname (string): '[output path]\[filename]'
%
% This code was taken from the Mathworks forum:
% https://www.mathworks.com/matlabcentral/answers/630163-how-can-i-create-
% a-3d-tiff-image-2d-stack-from-a-3d-matrix

for n = 1:size(stack, 4)
    if n == 1
        % First slice:
        imwrite(stack(:, :, :, n),fname)
    else
        % Subsequent slices:
        imwrite(stack(:, :, :, n),fname,'WriteMode','append');
    end 
end

end