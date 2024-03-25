function [wm_mask, gm_mask, tissue_mask] = create_wm_mask(mus,th)
%CREATE_WM_MASK create white matter mask from scattering map

%% Segment WM from mus map        
mus_th = mus;
mus_th(mus_th<=th) = 0;
mus_th = bwareaopen(mus_th, 50);

% Close the holes
se = strel('disk', 10);
close_mus = imclose(mus_th,se);

% Remove voxels with fewer than 20 connections
wm_mask = bwareaopen(close_mus, 20);

%%% Call the immask function within the MIMT toolbox
immask(wm_mask)

%% Segment GM from inverse of WM mask
% Create tissue mask from scattering mask, given the knowledge that the
% agarose scattering coeff. is approximately 0.
tissue_mask = mus;
tissue_mask(tissue_mask > 0.01) = 1;
se = strel('disk', 5);
tissue_mask = imerode(tissue_mask, se);
% Gray matter mask. Inverse of white matter mask and tissue mask
gm_mask = bsxfun(@times, tissue_mask, cast(~wm_mask,'like',tissue_mask));

%% Debugging figures of masks
figure; imagesc(mus_th); title('Thresholded mus');
figure; imagesc(close_mus); title('Thresholded/Closed');
figure; imagesc(wm_mask); title('WM Mask')
figure; imagesc(tissue_mask); title('Tissue')
figure; imagesc(gm_mask); title('GM Mask')


end