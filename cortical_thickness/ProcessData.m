% Process the OCT data
% Code by Iman Aganj.

%% Data folder
dataFolder = '/autofs/cluster/octdata2/users/Hui/ProcessHK001_cerebellum_20200618/Thickness/';
fileNames = {'seg_resample_mol.nii', 'seg_resample_gra.nii'};

%% Compile CortexMX
% mex MLI_thickness_MEX.c

% %% Make the line-integral masks
% L = makeLineSegments(60); %lSparse3(lSeg3DS(80, gridSphere(10), 5));
% save L_cerebellum L
% 
% % load('L_cerebellum');

%% Proceed with jobOCT
tr = .5;
nR = 1;
nTr = 0;
% n=1;
% parpool(length(fileNames))
% parfor n=1:length(fileNames)
for n=1:length(fileNames)
    disp([fileNames{n} ' ...'])
    jobOCT(dataFolder, fileNames{n}, tr, nR, nTr);
    
end

%% Combine the thickness and the skeleton
%r = 4;
% for n = 1:length(fileNames)
%     disp([fileNames{n} ':'])
%     V = load_nifti(fullfile(dataFolder, fileNames{n}));
%     Vs = load_nifti(fullfile(dataFolder, ['Skeleton_' fileNames{n}]));
%     Vt = load_nifti(fullfile(dataFolder, ['Thickness_' fileNames{n}]));
%     Vt = Vt.vol .* Vs.vol;
%     for r = 5:2:15
%         disp(['r = ' num2str(r) ' ...'])
%         Vs.vol = Vt .* (imerode(V.vol, strel(sph(r,r,r,r-1,2*r-1)), 'same')>.7);
%         save_nifti(Vs, fullfile(dataFolder, ['ThicknessOnSkeleton_r' num2str(r) '_' fileNames{n}]));
%     end
% end
% clear V Vs Vt Vts
