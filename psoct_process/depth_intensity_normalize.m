function [stack_norm] = depth_intensity_normalize(vol, iref)
%% Normalize the PS-OCT reflectance intensity across the image
% The intensity of the PS-OCT signal decreases in the z dimension. This
% leads to a heterogeneous contrast across the depth of each slice. The
% goal of depth intensity normalization is to correct for this decay.
%
% INPUTS
%   vol (double matrix): tissue OCT volume
%   iref (single array): subset of OCT volume to be used for noramlizing
%           the intensity
% OUTPUTS:
%   voln (single matrix): depth-intensity normalized volume

%% Calculate scaling factor from the intensity reference stack (iref)
% Use the first frame as the intensity reference
slice0 = iref(:,:,1);
intensity_0 = mean(slice0(:));

% Calculate mean of each frame in the reference stack
intensity_i = squeeze(mean(iref, [1,2]));

% Calculate the scaling factor to normalize the remaining frames to the
% same intensity as the initial frame.
r = intensity_i ./ intensity_0;

%% Normalize each slice
% Initialize matrix for storing normalized volume
stack_norm = zeros(size(vol));

% Convert the volume from a uint8 to a double array
vol = im2double(vol);

% Iterate over each layer in the volume
for ii = 1:size(vol,3)
    stack_norm(:,:,ii) = vol(:,:,ii) ./ r(ii);
end

% Scale the normalized volume from double back to uint8
stack_norm = rescale(stack_norm, 0, 255);
stack_norm = uint8(stack_norm);

end