%% This script is for testing the depth intensity normalization function
clear; clc;
% File paths to non-normalied volume
dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\Martinos_Datasets\caa_6\frontal\process_run2\mus\patches';
fname = 'patch_9.nii';
fullpath = fullfile(dpath, fname);

% Import NIFTI (.nii)
vol = load_nifti(fullpath);
vol = vol.vol;
figure; imshow(mat2gray(vol(:,:,1)));

%% Call normalization function
volnorm = depth_normalize(vol, 5);
volumeViewer(volnorm);

%% Save the output
% fout = fullfile(dpath, 'patch_9_norm.tif');
% segmat2tif(volnorm, fout);