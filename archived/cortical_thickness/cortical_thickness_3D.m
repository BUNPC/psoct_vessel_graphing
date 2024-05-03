%% Create the input image
% Simulate the soft-segmentation image 'I' as the space between two spheres.
I = TIFF2MAT('mask_infra.tif');
I(I(:)==255) = 1;
% interlopate in z
[Xq, Yq, Zq] = meshgrid(linspace(1,size(I,2),size(I,2)), linspace(1,size(I,1),size(I,1)), linspace(1,size(I,3),size(I,3)*15));
I_interp = interp3(single(I), Xq, Yq, Zq);

%% Create and save the line segments
% The radius parameter should be a few voxels larger than the expected
% maximum thickness. See the help of makeLineSegments.m for more parameters
% and parallelization options.
N = 85;
L = makeLineSegments(N);
save LineSegments L

%% Compute the thickness and the distance transform of the skeleton
% See the help of MLI_thickness.m for parameters and parallelization options.
% MLI_thickness_MEX.c must be compiled using: "mex MLI_thickness_MEX.c".
% Alternatively, the compiled files can be downloaded from:
% www.nitrc.org/projects/thickness
load LineSegments
[Thickness, SkeletonDistance] = MLI_thickness(I_interp, L);

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
