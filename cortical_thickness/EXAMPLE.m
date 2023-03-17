% This is an example script showing how to use this toolbox to compute the
% thickness and the skeleton of the tissue from a soft-segmented
% (probability) image, using minimum line integrals, as described in:
% 
% I. Aganj, G. Sapiro, N. Parikshak, S. K. Madsen, and P. Thompson,
% "Measurement of cortical thickness from MRI by minimum line integrals on
% soft-classified tissue," Human Brain Mapping, vol. 30, no. 10,
% pp. 3188-3199, 2009. http://doi.org/10.1002/hbm.20740
%
% See also:   MLI_thickness, makeLineSegments
%
% Codes by Iman Aganj.
% http://nmr.mgh.harvard.edu/~iman

%% Create the input image
% Simulate the soft-segmentation image 'I' as the space between two spheres.
N = 100;
[X,Y,Z] = ndgrid(-N/2:N/2);
I = (X.^2+Y.^2+Z.^2<(.4*N)^2) & ~(X.^2+(Y-.05*N).^2+Z.^2<(.3*N)^2);

%% Create and save the line segments
% The radius parameter should be a few voxels larger than the expected
% maximum thickness. See the help of makeLineSegments.m for more parameters
% and parallelization options.
L = makeLineSegments(round(.2*N));
save LineSegments L

%% Compute the thickness and the distance transform of the skeleton
% See the help of MLI_thickness.m for parameters and parallelization options.
% MLI_thickness_MEX.c must be compiled using: "mex MLI_thickness_MEX.c".
% Alternatively, the compiled files can be downloaded from:
% www.nitrc.org/projects/thickness
load LineSegments
[Thickness, SkeletonDistance] = MLI_thickness(I, L);

%% Visualize the results
% Volumetric visualization
figID = figure('units','normalized','outerposition',[0 0 1 1]); colormap jet
ThicknessOnSkeleton = Thickness.*(SkeletonDistance<1).*imerode(I, strel('sphere',2), 'same'); % (SkeletonDistance<1) is a 2-voxel-wide skeleton mask. Eroding 'I' eliminates border artifacts.
cmap = [0 max(Thickness(:))];
cmap_sd = [0 max(SkeletonDistance(:))];
for n = 1:size(I,3)
    set(gcf, 'Name', ['Frame #' num2str(n) ' out of ' num2str(size(I,3))])
    subplot(1,4,1), imagesc(I(:,:,n)), title Image
    caxis([0 1]), axis equal tight; colorbar
    subplot(1,4,2), imagesc(Thickness(:,:,n)), title Thickness
    caxis(cmap), axis equal tight; colorbar
    subplot(1,4,3), imagesc(SkeletonDistance(:,:,n)), title('Distance from Skeleton')
    caxis(cmap_sd), axis equal tight; colorbar
    subplot(1,4,4), imagesc(ThicknessOnSkeleton(:,:,n)), title('Thickness on the Skeleton')
    caxis(cmap), axis equal tight; colorbar
    drawnow
    %pause  % Uncomment to pause at each frame.
end
close(figID)

% Show the measured thickness on the outer surface of the skeleton.
figure, colormap jet
isosurface(ThicknessOnSkeleton>0, 0, Thickness)
patch(isosurface(I, 0), 'EdgeColor', 'none', 'FaceAlpha', .1);
camlight, lighting gouraud
axis equal tight; colorbar
title('Thickness on the outer surface of the skeleton')
