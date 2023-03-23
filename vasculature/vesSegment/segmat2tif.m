function segmat2tif(seg, fname)
%segmat2tif convert 3D matrix of segmentation to a 2D stack TIF file
% INPUTS:
%   seg (double matrix): 3D matrix of vessel mask
%   fname (string): '[output path]\[filename]'
%
% This code was taken from the Mathworks forum:
% https://www.mathworks.com/matlabcentral/answers/630163-how-can-i-create-
% a-3d-tiff-image-2d-stack-from-a-3d-matrix


for n = 1:size(seg, 3)   
    if n == 1
        % First slice:
        imwrite(seg(:, :, n),fname)
    else
        % Subsequent slices:
        imwrite(seg(:, :, n),fname,'WriteMode','append');
    end 
end

end