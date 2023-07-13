%% Main script to overlay segmentation and original volume.

%% Initialize data path for linux or personal machine
% Check if running on local machine for debugging or on SCC for processing
if ispc
    %%% Local machine
    % Top-level directory
    dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\';
    % Directory to ps-oct volume
    voldir = '\dist_corrected\volume\';
    % Directory to segmentation
    segdir = '\dist_corrected\volume\gsigma_1-2-3-4-5_gsize_5--9-13-17-21\';
    % Subject IDs
    subid = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424',...
             'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_8572', 'CTE_7126'};
    % Filename to parse (this is test data)
    fname = 'ref_4ds_norm_inv_segment_pmin';
    % Probability matrix range
    pmin = 0.20 : 0.01 : 0.26;
    % filename extension
    ext = '.tif';
elseif isunix
    %%% Computing cluster (SCC)
    % Path to top-level directory
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_10T/';

end

%% Overlay
%%% Load PSOCT volume (TIF) and convert to MAT
for ii = 1:length(subid)
    % Define entire filepath 
    fullpath = fullfile(dpath, subid{ii}, segdir);
    % Define filename of original ps-oct volume
    filename = strcat(fullpath, strcat(fname, ext));
    % Convert .tif to .MAT
    vol = TIFF2MAT(filename);
    
    %%% Load segmentation across range of probabilities
    for j=1:6
        
    end

end




