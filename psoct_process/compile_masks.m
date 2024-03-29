%% Compile the masks
% Ayman created masks of the white matter and gray matter from the
% scattering coefficient maps. Each subject has a single mask for each
% scattering map. Since a tissue volume mask already exists for each
% subject, the simplest path forward is to perform a bitwise operation
% between the white matter mask and the tissue mask. This will also enable
% the generation of a gray matter mask. The masks are stored in:
% /projectnb/npbssmic/ns/Ann_Mckee_samples_55T/ROIs_ayman/result/...
%   David_result_ayman_GM_WM/Masks/
%
% The following work must be performed to apply these masks to the OCT
% intensity volumes:
% - replicate each mask 11 times for each physical slice (11 depths/slice)
% - combine the masks in sequential order
% - bitwise multiply the white matter mask with the tissue mask
% - bitwise multiply the inverse of the WM mask with the tissue mask. This
%       will generate the gray matter mask
%
% These subjects have masks:
% AD_10382, AD_20832, AD_20969, AD_21354, AD_21424,
% CTE_6489, CTE_6912, CTE_7019, CTE_7126,
% NC_6047, NC_6839, NC_7597, NC_8095

%% Initialize Environment
clear; clc; close all;

%%% Set maximum number of threads
% Retrieve the number of available cores
n_cores = str2double(getenv('NSLOTS'));
% Set the maximum number of threads equal to the number of cores
maxNumCompThreads(n_cores);

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

%% Initialize directories and filenames
%%% Directories and file name
% Upper level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T';    
% Subfolder containing ref files
subdir = '/dist_corrected/volume/ref';
% non-normalized PSOCT volume
vol_name = 'ref';
% Directory containing masks
mdir = fullfile(dpath,'/ROIs_ayman/result/David_result_ayman_GM_WM/Masks/');

%% Compile masks

%%% Subjects with masks
subid = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_7126',...
         'NC_6047', 'NC_6839', 'NC_7597', 'NC_8095'};
% Number of physical slices for each subject
nslice = 10;


%%% Iterate through subjects

for ii = 1:length(subid)
    % Load the OCT ref file
    
    % Initialize the matrix for storing the mask
    
    %%% Iterate over each slice and create a mask
    
    % Load the mask from Ayman
    
    % Replicate the mask

    % Add mask to the matrix


end





































