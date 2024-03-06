function [stack_norm] = depth_intensity_normalize(vol, kernel1)
%% Normalize the PS-OCT reflectance intensity across the image
% The intensity of the PS-OCT signal decreases in the z dimension. This
% leads to a heterogeneous contrast across the depth of each slice. The
% goal of depth intensity normalization is to correct for this decay.

% This code was originally written by Stephan Zhang and then adapted by
% Mack Hyman. Stephan's code convolves one slice of the image to generate
% the local average intensity. Then, he normalizes the slice by this
% convolved slice. However, one issue is that this just normalizes each
% slice to itself rather than to the intensity decay. The results still
% contain a trace line at the interface between slices, so the original
% issue was not fully resolved.

% Further work:
%   - normalize across the z dimension for each layer of each physical slice.
%   - smooth the boundary between slices

%% Initialization
% Convert volume to single
vol = single(vol);
stack_norm = zeros(size(vol));

%% Normalize each slice

% Iterate over each layer in the volume
for ii = 1:size(vol,3)
    slice = vol(:,:,ii);
%         slice(slice<3000)=0;
%         tmp=mean(slice(slice>0.1));
    tmp = single(convn(slice, ones(kernel1,kernel1)./kernel1^2,'same')+500);
    stack_norm(:,:,ii)=single(slice./tmp);
end

%% Convert stack to unsigned int 16
stack_norm = uint16(stack_norm./2*65535);

end