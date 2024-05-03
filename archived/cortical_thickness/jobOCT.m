% Compute the thickness and the skeleton for one OCT probability image
% Code by Iman Aganj.

function jobOCT(dataFolder, fileName, tr, nR, nTr)

load L_cerebellum

V = load_nifti(fullfile(dataFolder, fileName));
I = V.vol;

parpool(6)
param.threshStop = nR;
param.numVoxStop = tr;
param.numVoxValey = nTr;
param.useParpool = true;
[Thickness, SkeletonDistance] = MLI_thickness(I, L, param);
% [J, ~, Ja] = CortexMX(V.vol, L, nR, tr, nTr);
clear L

V.vol = SkeletonDistance;
save_nifti(V, fullfile(dataFolder, ['SkeletonDistance_' fileName]));

V.vol = Thickness;
save_nifti(V, fullfile(dataFolder, ['Thickness_' fileName]));

% im_filter = imerode(V.vol, strel('ball',10,10), 'same');
% ThicknessOnSkeleton = Thickness.*(SkeletonDistance<1).*(im_filter>(min(im_filter(:)+0.7))); % (SkeletonDistance<1) is a 2-voxel-wide skeleton mask. Eroding 'I' eliminates border artifacts.
% V.vol = ThicknessOnSkeleton;
% save_nifti(V, fullfile(dataFolder, ['ThicknessOnSkeleton_' fileName]));
clear Thickness SkeletonDistance ThicknessOnSkeleton
poolobj = gcp('nocreate');
delete(poolobj);


