%% Main script to create the white matter mask
% Each tissue volume has disparate sucli and gyri masks. They are rough
% croppings of the sulci and gyri, so they need to be combined with the
% tissue masks to refine their borders. This script will combine the
% disparate sucli masks into a single sulci mask. It will do the same for
% the gyri. Then, it will combine the sulci and gyri with the tissue mask.

% Remaining work:
%{
- combine the sulci into single mask
- combine the gyri into single mask
- Mask the sulci and gyri (with tissue mask)
%}
clear; clc; close all;

%% Add top-level directory of code repository to path
% Print current working directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Remove the two sub folders to reach parent
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));

%% Set maximum number of threads to avoid session termination
% Retrieve the number of available cores
n_cores = str2double(getenv('NSLOTS'));
% Set the maximum number of threads equal to the number of cores
maxNumCompThreads(n_cores);

%% Initialize directories and filenames
%%% Directories and file name
% Upper level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T';    
% Subfolder containing ref files
subdir = '/dist_corrected/volume/ref';
% non-normalized PSOCT volume
vol_name = 'ref1';

%% Debugging Variables
% Debugging figures visualization
viz = false;
% Boolean to load in the cleaned tissue masks
tmask_bool = true;

%% Initialize subject ID lists
%%% All subjects to analyze
subid = {'AD_10382', 'AD_20832', 'AD_20969',...
         'AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912',...
         'CTE_7019', 'CTE_7126','CTE_8572',...
         'NC_6839',  'NC_6974', 'NC_8653',...
         'NC_21499', 'NC_301181'};

%%% stacks that require truncating
% Last depth to retain for each stack
zmins = [187, 165, 242, 165, 220,...
        242, 110, 198, 198, 198,...
        176, 198, 165, 165, 220];
% Create dictionary to store last image in stack
d = dictionary(subid, zmins);

%% Create white mask for subjects
for ii = 1:length(subid)
    %%% Debugging information
    fprintf('\n---------Starting Subject %s---------\n',subid{ii})

    %%% Import tissue mask to get dimensions
    fprintf('importing non-normalized volume\n')
    fpath = fullfile(dpath, subid{ii}, subdir, '/masks/','mask_tiss.mat');
    mask_tiss = load(fpath); 
    mask_tiss = mask_tiss.mask_tiss;

    %%% Retrieve last slice # and N slices per physical stack
    sid = char(subid{ii});
    zmax = d(subid(ii));
    
    %%% Initialize matrices for storing masks
    mask_sulci = false(size(mask_tiss,1),size(mask_tiss,2),zmax);
    mask_gyri = false(size(mask_tiss,1),size(mask_tiss,2),zmax);
    
    %%% Initialize directory for storing masks
    mdir = fullfile(dpath,sid,subdir,'/masks/');

    %% Combine sulci TIFs into single sulci mask
    % Extract file names of all sulcus TIFs in the mask folder
    files = dir(fullfile(mdir, 'sulcus*'));
    files = {files.name}.';
    % Load each and combine into single matrix
    for j = 1:length(files)
        % Load in sulcus mask TIF
        tmp = TIFF2MAT(fullfile(mdir,files{j}));
        tmp = logical(tmp);
        % Add to main sulcus mask matrix
        mask_sulci = mask_sulci + tmp(:,:,1:zmax);
    end

    %% Combine gyri TIFs into single sulci mask
    % Extract file names of all sulcus TIFs in the mask folder
    files = dir(fullfile(mdir, 'gyrus*'));
    files = {files.name}.';
    % Load each and combine into single matrix
    for j = 1:length(files)
        % Load in sulcus mask TIF
        tmp = TIFF2MAT(fullfile(mdir,files{j}));
        tmp = logical(tmp);
        % Add to main sulcus mask matrix
        mask_gyri = mask_gyri + tmp(:,:,1:zmax);
    end

    %% Combine sulci/gyri masks with tissue mask
    % Load in the tissue mask
    mask_tiss = load(fullfile(mdir,"mask_tiss.mat"));
    mask_tiss = mask_tiss.mask_tiss;
    mask_tiss = mask_tiss(:,:,1:zmax);

    % Bitwise multiply tissue mask with sulci and gyri masks
    mask_sulci = bsxfun(@times, mask_tiss ,cast(mask_sulci,'like',mask_tiss));
    mask_sulci = logical(mask_sulci);
    mask_gyri = bsxfun(@times, mask_tiss ,cast(mask_gyri,'like',mask_tiss));
    mask_gyri = logical(mask_gyri);

    %% Save Masks
    % Save outputs as .MAT
    fprintf('Saving outputs as .MAT\n')
    sulci_output = fullfile(mdir,'mask_sulci.mat');
    gyri_output = fullfile(mdir,'mask_gyri.mat');
    save(sulci_output,"mask_sulci",'-v7.3');
    save(gyri_output,"mask_gyri",'-v7.3');
    
    % Save outputs as .TIF
    fprintf('Saving outputs as .TIF\n')
    sulci_output = fullfile(mdir,'mask_sulci.tif');
    gyri_output = fullfile(mdir,'mask_gyri.tif');
    segmat2tif(im2uint8(mask_sulci),sulci_output)
    segmat2tif(im2uint8(mask_gyri),gyri_output);
    
    % Save outputs as .NII
    res = [0.012, 0.012, 0.015];
    dtype = 'uchar';
    permuteflag = 1;
    fprintf('Saving outputs as .NII\n')
    save_mri(mask_sulci, fullfile(mdir,'mask_sulci.nii'), res, dtype, permuteflag);
    save_mri(mask_gyri, fullfile(mdir,'mask_gyri.nii'), res, dtype, permuteflag);

end