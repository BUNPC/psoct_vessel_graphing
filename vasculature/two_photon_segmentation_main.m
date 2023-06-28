%% Main script for segmenting 2PM images
% Author: Mack Hyman
% Email: mhyman@bu.edu
% Date Created: May 15, 2023
%
% Detailed Description
%{
This script performs the following:
- segment the original volume
- apply a mask to the segmentation
%}
clear; clc; close all;

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

%% Import volume (.TIF or .BTF) & convert to MAT 

% Check if running on local machine for debugging or on SCC for processing
if ispc
    %%% Local machine
    dpath = ['C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\'...
        'test_data\BA4445_samples_16T'];
    % Subject IDs
    subid = {'\BA4445_I38_2PM'};
    subdir = '\aip\';
    % Number of slices to parse
    slices = [10];
    % Filename to parse (this is test data)
    fname = 'channel1';
    % filename extension
    ext = '.tif';
elseif isunix
    %%% Computing cluster (SCC)
    % Path to top-level directory
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_10T/';
    % Subject IDs
%     subid = {'AD_10382_2P', 'AD_20832_2P', 'AD_20969_2P', 'AD_21354_2P', 'AD_21424_2P',...
%              'CTE_6489_2P', 'CTE_6912_2P', 'CTE_7019_2P', 'CTE_8572_2P', 'CTE_7126_2P',...
%              'NC_21499_2P', 'NC_6047_2P', 'NC_6839_2P', 'NC_6974_2P', 'NC_7597_2P',...
%              'NC_8095_2P', 'NC_8653_2P'};
    subid = {'AD_10382_2P', 'AD_20832_2P'};
    subdir = '/aip/';
    % Filename to parse (this will be the same for each subject)
    fname = 'channel1';
    % filename extension
    ext = '.tif';
end

%% Initializatoin
%%% Post-processing parameters
% Minimum number of voxels to classify as segment
% A segment with < "min_conn" voxels will be removed
min_conn = 30;

% Length of line in the structuring element (strel) for expanding lines
line_len = 3;

% Disk radius
dr = 3;

%%% Preliminary file paths
meth = {'Sobel', 'Prewitt', 'Roberts', 'log', 'zerocross', 'Canny'};
% Define filename for slice
slice_fname = 'channel1-10_16b_cropped';
% Define entire filepath 
subpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\BA4445_samples_16T\BA4445_I38_2PM\aip\';
% Filepath for input .TIF
input = fullfile(subpath, strcat(slice_fname,'.tif'));
% Convert .tif to .MAT
slice = TIFF2MAT(input);
% sub-directory for filtered outputs
slicepath = fullfile(subpath, strcat('slice_10'));

%%% FFT
slice_freq = fftshift(fft2(slice));
freq_amp = log(abs(slice_freq));

%% Frequency filters (LPF, HPF)
%%% Plot FFT prior to filtering
figure; imshow(freq_amp, [])
figure;
subplot(3,1,1); imshow(freq_amp, []); axis on; title('FFT of slice')

%%% Create circular mask (LPF and HPF)
rows = size(slice_freq, 1);
cols = size(slice_freq, 2);
[x, y] = meshgrid(1:cols, 1:rows);
% Create matrix of zeros
lpf = zeros(rows, cols);
c_x0 = cols./2;
c_y0 = rows./2;
r = 500;
% Define inside circle as 1 (low pass filter)
lpf((x - c_x0).^2 + (y - c_y0).^2 <= r.^2) = 1; 
% Define high pass filter
hpf = ~lpf;

%%% Apply LPF & HPF filters
% Multiply LPF frequency filter with FFT of slice
slice_freq_lpf = lpf .* slice_freq;
subplot(3,1,2); imshow(log(abs(slice_freq_lpf)),[]); axis on; title('LPF')
% Multiply HPF frequency filter with FFT of slice
slice_freq_hpf = hpf .* slice_freq;
subplot(3,1,3); imshow(log(abs(slice_freq_hpf)),[]); axis on; title('HPF')

%%% Plot original and filtered images
figure;
subplot(3,1,1); imshow(slice, []); axis on; title('Original Slice')
% Convert filtered frequency images back to spatial domain
slice_lpf = real(ifft2(ifftshift(slice_freq_lpf)));
slice_hpf = real(ifft2(ifftshift(slice_freq_hpf)));
subplot(3,1,2); imshow(slice_lpf, []); title('LPF image')
subplot(3,1,3); imshow(slice_hpf, []); title('HPF image')

%% Filter with elipsoid

%% Filter noise peaks above threshold
% Find frequency components w/ amp > threshold
th = 16.5;
fpeaks = freq_amp > th;
figure;
subplot(1,2,1); imshow(freq_amp, []); title('FFT')
subplot(1,2,2); imshow(fpeaks, []); title('Peaks of FFT')

%%% Add rectangle to mask
% Initialize ones matrix size of freq matrix
noise_mask = ones(size(slice_freq));
% Add zeros at higher intensities at high frequencies
% Left side
noise_mask(2838:3004, 1:188) = 0;
noise_mask(4724:4900, 1:650) = 0;
noise_mask(3796:3954, 634:931) = 0;
% Right side
noise_mask(2852:3000, 10039:10608) = 0;
noise_mask(4724:4863, 10420:10608) = 0;
noise_mask(3796:3954, 9666:10006) = 0;

figure; imagesc(noise_mask); axis on;

%%% Apply filter to FFT
slice_filt = noise_mask .* slice_freq;
figure;
slice_filt_real = real(ifft2(ifftshift(slice_filt)));
imshow(slice_filt_real,[])



%% Iterate through edge detection filters
for ii = 1:length(meth)
    % Segment
    volmask = edge(slice, meth{ii});
    % Clean + dilate segmentation
    volmask_clean = clean_dilate(volmask, min_conn, dr, subpath, slice_fname);

    %%% Save outputs (segmented & cleaned)
    % Base name
    proc = strcat(slice_fname,'_',meth{ii});
    % save segmented
    filename = fullfile(slicepath, strcat(proc, ext));
    imwrite(uint16(volmask), filename);
end

%% Frangi 2P Segmentation (Mathworks fibermetric)

%%% Frangi Filter Parameters
% vessel thickness (pixels). default: [4 6 8 10 12 14]
thick = 8;
% StructureSensitivity: vessel contrast threshold.
% default = 0.01*diff(getrangefromclass(I))
sense = 2;
% ObjectPolarity: 'bright' (default) or 'dark'
op = {'bright', 'dark'};

%%% Iterate subjects
for ii = 1:length(subid)
    %%% Iterate slices
    for j = 1:length(slices)
        % Define filename for slice
        slice_fname = strcat(fname,'-',num2str(slices(j)),'_16b_cropped');
        % Define entire filepath 
        subpath = fullfile(dpath, subid{ii}, subdir);
        filename = strcat(subpath, strcat(slice_fname, ext));
        % Make sub-directory for slice
        slicepath = fullfile(subpath, strcat('slice_', num2str(slices(j))));
        if not(isfolder(slicepath))
            mkdir(slicepath)
        end
        % Convert .tif to .MAT
        slice = TIFF2MAT(filename);
        for k = 1:length(thick)
            for lpf = 1:2
                %%% Frangi filter + clean
                % Segment
                volmask = fibermetric(slice, thick(k),'ObjectPolarity',op{lpf});
                % Clean + dilate segmentation
                volmask_clean = clean_dilate(volmask, min_conn, dr, subpath, slice_fname);
        
                %%% Save outputs (segmented & cleaned)
                % Base name
                proc = strcat(slice_fname,'_frangi_',op{lpf},'_',num2str(thick(k)));

                % save segmented
                filename = fullfile(slicepath, strcat(proc, ext));
                imwrite(uint16(volmask), filename);
        
                % save segmented + cleaned
                tmp_name = strcat(proc,'_cleaned');
                filename = fullfile(slicepath, strcat(tmp_name, ext));
                imwrite(uint16(volmask_clean), filename);
            end
        end
    end
end


%% Frangi 2P Segmentation (custom scripts)
%{
%%% 2P microscopy voxel will always be [5, 5] micron
vox_dim = [2,2];

%%% Std. Dev. for gaussian filter (range)
% The value of sigma corresponds to the smallest resolvable radius of the
% vessels. If sigma equals 1, then the smallest resolvable vessel will be
% 1*voxel. For 2PM, the smallest resolvable vessel has a radius = 5um
sigma = [1, 3];

%%% Minimum number of voxels to classify as segment
% A segment with < "min_conn" voxels will be removed
min_conn = 30;

%%% Lenght of line in the structuring element (strel) for expanding lines
line_len = 3;

%%% Boolean for converting segment to graph (0 = don't convert, 1 = convert)
graph_boolean = 0;

%%% Segment each slice for each subject
for ii = 1:length(subid)
    for j = 1:length(slices)
        % Define filename for slice
        slice_fname = strcat(fname,'-',num2str(slices(j)),'_cropped_inv');
        % Define entire filepath 
        subpath = fullfile(dpath, subid{ii}, subdir);
        filename = strcat(subpath, strcat(slice_fname, ext));

        % Make sub-directory for slice
        slicepath = fullfile(subpath, strcat('slice_', num2str(slices(j))));
        if not(isfolder(slicepath))
            mkdir(slicepath)
        end
        % Convert .tif to .MAT
        slice = TIFF2MAT(filename);
        % Perform segmentation
        [seg, fname_seg] = segment_main(slice, sigma, slicepath, slice_fname);

        %%% Filter output of Frangi filter
        for k=1:length(line_len)
            % Remove small segments and apply mask
            [seg_masked, fout] = ...
                clean_dilate(seg, min_conn, line_len, slicepath, slice_fname);
            %%% Convert masked segmentation to graph
            if graph_boolean
                % Use masked segmentation to create graph
                Graph = seg_to_graph(seg_masked, vox_dim);
                
                % Initialize graph information (work in progress)
                % Graph = initialize_graph(Graph);
        
                % Create new filename for graph and add .MAT extension
                fname_graph = strcat(fname_seg,'_mask', num2str(radii(k)),'_graph.mat');
                fout = strcat(slicepath, fname_graph);
                save(fout,'Graph');
            end
        end
    end
end
%}

%% Apply Mask
function [seg_final] =...
    clean_dilate(seg, min_conn, dr, fullpath, fname, vth)
% Remove the edges labeled as vessels.
%   INPUTS:
%       seg (matrix) - output of segmentation function
%       min_conn (double) - minimum number of connections to classify as
%                           vessel
%       line_len (double) - length of line for strel
%       fullpath (string) - absolute directory for saving processed data
%       fname (string) - filename prior to applying mask
%       vth (int) - width threshold in pixels
%   OUTPUTS:
%       seg_masked (matrix) - seg with boundaries eroded to remove
%           erroneously labeled vessels.

%%% Remove segments with fewer than "min_conn" connected voxels
figure; imshow(seg); title('Frangi - before processing')
[seg, ~, ~] = rm_small_seg(seg, min_conn);
figure; imshow(seg); title('Frangi - small segments removed')

%%% Dilate the lines
%%% Dilate with perpendicular lines
% Create perpendicular structure
% se90 = strel('line',line_len,90);
% se0 = strel('line',line_len,0);
% Dilate gradient mask
% seg_dil = imdilate(seg, [se90 se0]);

%%% Dilate with disk
se = strel("disk", dr);
seg_dil = imdilate(seg, se);
figure; imshow(seg_dil); title('Frangi - Dilated Gradient Mask')

%%% Fill interior gaps
seg_fill = imfill(seg_dil,'holes');
figure; imshow(seg_fill); title('Frangi - Filled Holes')

%%% Smooth segmentation
seD = strel('diamond',1);
seg_final = imerode(seg_fill,seD);
figure; imshow(seg_fill); title('Frangi - Eroded')

% figure; 
% bw = imbinarize(seg_final);
% figure; imshow(labeloverlay(seg_final, bw)); title('Frangi - Final')

%%% Remove small segments
[seg_final, ~, ~] = rm_small_seg(seg_final, min_conn);
figure; imshow(seg_fill); title('Frangi - Small Segments Removed')

%%% Save segmented/masked volume as .MAT and .TIF
%{
% Convert processed matrix to tif
tmp_fname = strcat(fname,'_frangi',vth);
fout = strcat(fullpath, tmp_fname, '.tif');
segmat2tif(seg_cleaned, fout);

% Save vessel segment stack as .MAT for the next step (graphing)
fout = strcat(fullpath, tmp_fname, '.mat');
save(fout, 'seg_cleaned', '-v7.3');
%}
end

%% Initialization of vesGraphValidate
function [graph_init] = initialize_graph(Graph)
%%% Perform the manual operations for initializing data in the GUI.
% Run "Verification > get segment info > Update"
% Run "Update branch info"
% Run "Regraph Nodes" to down sample
% Open GUI with both image and data (graph)
% Run prune_loops and prune_segment
% Run straighten
graph_init = Graph;
end

%% Segment volume
function [seg, fname] =...
    segment_main(slice, sigma, fullpath, fname)
% Multiscale vessel segmentation
%   INPUTS:
%       slice (matrix) - the original volume prior to segmentation
%       FrangiScaleRange (array) - range of std. dev. values of gaussian
%            filter to calcualte hessian matrix at each voxel
%       min_prob - threshold to determine which voxel belongs to a vessel.
%           Applied to probability matrix from frangi filter output
%       min_conn - vesSegment uses the function bwconncomp to determine the
%               number of connected voxels for each segment. If the number
%               of voxels is less than this threshold, then the segment
%               will be removed.
%       fullpath (string) - absolute directory for saving processed data
%       fname (string) - filename of segmentation
%       ext (string) - filename extension (.TIF or .MAT)
%   OUTPUTS:
%       seg - segmentation of vessels (Frangi filter of original volume)
%       fname (string) - upadted filename prior to applying mask

%% convert volume to double matrix
I = double(slice);

%%% Range of sigma values
options.FrangiScaleRange = sigma;

%%% Segment the original volume
[seg, ~, ~] = FrangiFilter2D(I, options);

%%% Save segmentation
fname = strcat(fname,'_segment','_sigma_',...
    num2str(sigma(1)), '_', num2str(sigma(2)));
fout = fullfile(fullpath, strcat(fname, '.tif'));
segmat2tif(seg, fout);

end