function[voln]=depth_normalize_2(vol, kernel, th)
% DEPTH_NORMALIZE_2 normalize the depth intensity attenuation.
% This was written by Stephan Chang.
% INPUTS:
%   vol (uint16 matrix): non-normalized OCT volume. Either the "ref_10ds.btf" or
%               "ref_4ds.btf"
%   kernel (uint8): one dimensiona of a 2D kernel for normalization. This
%                   value is typically set to 50. It is later used to
%                   create a blurring matrix of ones
%                   (eg ones(kernel, kernel).
% OUTPUTS:
%   voln (matrix): normalzied volume

%% Depth Intensity Normalization

% Initialize matrix for storing normalized volume
voln = uint16(zeros(size(vol)));

% Initialize a blurring kernel for convolution
kernel = single(ones(kernel,kernel))./(single(kernel).^2);

% Iterate over each depth in the volume
for ii=1:size(vol,3)
    % Take the ii slice
    depth = single(vol(:,:,ii));
%     figure; imagesc(depth); colorbar; title('depth before threshold');
    % apply minimum threshold
    depth(depth<th) = 0;
%     figure; imagesc(depth); colorbar; title('depth after threshold');

    % Create normalization matrix (convolve depth with kernel, add offset)
    norm_mat = single(convn(depth, kernel, 'same') + 500);
    
    % Store normalized depth in matrix
    voln(:,:,ii) = im2uint16(depth./norm_mat);
%     close all;
end

end