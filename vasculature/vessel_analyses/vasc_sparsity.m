function [vs] = vasc_sparsity(seg, mask, vox_dim)
%VASC_SPARSITY compute the distance of each voxel from a blood vessel
%vessel
%   Import the vessel segmentation, the masked OCT volume, and the voxel
%   dimensions. For each point in the masked OCT volume, find the minimum
%   distance to the nearest blood vessel. 
%
%   INPUTS:
%       seg (logical matrix): segmentation matrix
%       mask (logical matrix): tissue mask for region
%       vox_dim (double array): dimensions of 
%   OUTPUTS:
%       vs (double array): Euclidean distance of each non-vessel voxel from
%           the nearest blood vessel (microns)
%

%% Identify indices of non-vessel voxels in the masked OCT matrix
% non-zero indices in the mask
pause(0.10);

% non-zero indices in segmentation

%% Find distance to nearest vessel


%% Output results

end