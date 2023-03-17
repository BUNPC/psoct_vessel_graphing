%% batch job version of vessel segmentation
% Author: Jiarui Yang
% 10/21/20
% Walkthrough:
%{
To Run this script, open a terminal in SCC and navigate to:
/projectnb/npbssmic/s/Matlab_code/vesSegment

Run the command: 
qsub vesSeg.sh
%}

% To Do:
%{
- Add relative path
- Determine input parameters for vesSegment
- 

%}

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
% (psoct_human_brain\vasculature\vesSegment)
newdir = mydir(1:idcs(end-1));
addpath(genpath(newdir));

%% load volume
dpath = '/projectnb/npbssmic/ns/Hui_Frangi_dataset/200218depthnorm/';
fname = 'volume_ori_inv';
filename = strcat(dpath, strcat(fname, '.tif'));
ref = TIFF2MAT(filename);

%% multiscale vessel segmentation
% I = 3D angiogram
% sigma = vector of standard deviation values of gaussian filter to
% thres = threshold value to get segmentation [0, 1]
I = double(ref);
sigma = 1;
thres = 0.2;
[~, I_seg]=vesSegment(I, sigma, thres);

%% save segmentation
% Save vessel segment stack as .MAT for the next step (graph recon)
fout = strcat(dpath, fname, '_sigma', num2str(sigma));
save(strcat(fout, '.mat'), 'I_seg', '-v7.3');

%% Save as .TIF for visualizations
tifout = strcat(fout, '.tif');
segmat2tif(I_seg, tifout);
