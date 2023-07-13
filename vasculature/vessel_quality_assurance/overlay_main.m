%% Main script to overlay segmentation and original volume.
%{
The purpose of this script is to overlay the vessel segmentation and the
PS-OCT volume (each slice).
%}

%% Add top-level directory of code repository to path
% This allows Matlab to find the functions in the project folders

% Start in current directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Truncate path to reach top-level directory (psoct_vessel_graphing)
topdir = mydir(1:idcs(end));
addpath(genpath(topdir));

%% Initialize data path for linux or personal machine
% Check if running on local machine for debugging or on SCC for processing
if ispc
    % Top-level directory
    dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\';
    
    % Volume directory + volume filename (same for each subject)
    voldir = '\dist_corrected\volume\';
    volname = 'ref_4ds_norm_inv';

    % Directory to segmentation
    segdir = '\dist_corrected\volume\gsigma_1-2-3-4-5_gsize_5--9-13-17-21\';
    % Subject IDs
    subid = {'NC_6839'};
    
    % Filename to parse (this is test data)
    segname = 'ref_4ds_norm_inv_segment_pmin_0.23_mask40';
    % Probability matrix range
    pmin = 0.20 : 0.01 : 0.26;
    % filename extension
    ext = '.tif';

elseif isunix
    % Path to top-level directory
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_10T/';

end

%% Overlay
%%% Load PSOCT volume (TIF) and convert to MAT
for ii = 1:length(subid)
    %%% Load segmentation
    % Define entire filepath 
    fullpath = fullfile(dpath, subid{ii}, segdir);
    % Define filename of original ps-oct volume
    filename = strcat(fullpath, strcat(segname, ext));
    % Convert .tif to .MAT
    seg = TIFF2MAT(filename);
    
    %%% Load PS-OCT volume
    fullpath = fullfile(dpath, subid{ii}, voldir);
    % Define filename of original ps-oct volume
    filename = strcat(fullpath, strcat(volname, ext));
    % Convert .tif to .MAT
    vol = TIFF2MAT(filename);
    pause(0.1)

    %%% Overlay segmentation and probaiblity
    tmp = imoverlay(vol(:,:,1), seg(:,:,1), 'green');
    ov = overlay_vol_seg(vol, seg, 'green');
    % Save output volume
    fullpath = fullfile(dpath, subid{ii}, segdir);
    % Define filename of original ps-oct volume
    filename = strcat(fullpath, strcat(volname, '_overlay', ext));
    segmat2tif(ov, filename)
    
    %%% Overlay segmentation across range of probabilities

end




