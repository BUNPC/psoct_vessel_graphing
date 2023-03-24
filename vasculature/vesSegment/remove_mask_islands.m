function [mask_rm] = remove_mask_islands(mask, range)
%REMOVE_MASK_ISLANDS Remove small groups from mask matrix
% The mask contains some small "islands" of nonzero elements. These are
% erroneously labeled as nonzero and need to be removed.
%%% INPUTS:
%       mask (3D uint16): zero for areas to remove and nonzero to keep.
%       range (uint16): the function bwareafilt finds groups within a 2D
%                       matrix and then calculates the size of the object
%                       in pixels. Groups within this range will be kept
%                       and groups smaller or larger will be discarded.
%                       Here range = [min, max]
%%% OUTPUTS:
%       mask_rm (3D uint16): the original mask will islands removed

% Reshape 3D matrix into 2D matrix
[mask_re, r, c, s] = reshape_3dmat_to_2d(mask);

% Convert to logical to pass into bwareafilt
mask_re = logical(mask_re);

% Find/delete islands in the mask
mask_rm = bwareafilt(mask_re, range);

% Debugging figure
figure;
subplot(2,1,1); imshow(mask_re(:,1:(2*c))); title('Original Mask')
subplot(2,1,2); imshow(mask_rm(:,1:(2*c))); title('Cleaned Mask')


% Convert back to 3D matrix
mask_rm = reshape(mask_rm, [r, c, s]);

end

