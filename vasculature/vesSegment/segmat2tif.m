function segmat2tif(seg, fname)
%segmat2tif convert 3D matrix of segmentation to a 2D stack TIF file
% INPUTS:
%   seg (double matrix): 3D matrix of vessel mask
%   fname (string): '[output path]\[filename]'
%
% This code was taken from the Mathworks forum:
% https://www.mathworks.com/matlabcentral/answers/630163-how-can-i-create-
% a-3d-tiff-image-2d-stack-from-a-3d-matrix


% Normalize for RGB colormap mapping:
% scaled_seg = seg./max(seg, [], 'all');

% Colormap
% cmp = parula;

% Write flag: run once with iwrite = false if you want to see what you are writing first
iwrite = true;

for n = 1:size(seg, 3)
    
    % Make an RGB image:
%     i_img = ind2rgb(round(scaled_seg(:, :, n)*256), cmp);
    
    % Check what you are writing
%     figure; imagesc(seg(:, :, n));
    
    % Generate your tiff stack:
    if iwrite
        if n == 1
            % First slice:
            imwrite(seg(:, :, n),fname)
        else
            % Subsequent slices:
            imwrite(seg(:, :, n),fname,'WriteMode','append');
        end 
    end
    
    disp(n)
end

end