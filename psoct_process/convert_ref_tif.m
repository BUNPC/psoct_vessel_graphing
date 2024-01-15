%% Combine the b-scans (ref#.mat) to remake the c-scan
% Each subdirectory ".../[subjectID]/dist_corrected/volume/ contains the
% b-scans in the format Ref_BASIC*#.btf. These need to be converted to .MAT
% format so that they can be compiled into a c-scan and used to generate a
% tissue and pia mask.
% Subjects and respective filenames that need to be converted:
% Subject ID: NC_301181 -- Ref_BASIC_4ds#.btf
% Subject ID: NC_7597 -- Ref_BASIC#.btf
% Subject ID: NC_8095 -- Ref_BASIC#.btf

% TODO:
% - verify correct file paths
% - ensure "ref.mat" and "ref.tif" not in output folders

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

%% Initialize data path for linux or personal machine (debugging)

%%% Computing cluster (SCC)
if isunix
    % Path to top-level directory
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
    % Set # threads = # cores for job
%     NSLOTS = str2num(getenv('NSLOTS'));
%     maxNumCompThreads(NSLOTS);
% 
%     % Check to see if we already have a parpool, if not create one with
%     % our desired parameters
%     poolobj = gcp('nocreate');
%     if isempty(poolobj)
% 	    myCluster=parcluster('local');
% 	    % Ensure multiple parpool jobs don't overwrite other's temp files
% 	    myCluster.JobStorageLocation = getenv('TMPDIR');
% 	    poolobj = parpool(myCluster, NSLOTS);
%     end
else
    error('Should only be run on UNIX')
end

% Subfolder containing data
subdir = '/dist_corrected/volume/';
% Subject IDs
subid = {'NC_7597','NC_8095','NC_301181'};


for ii = 1:length(subid)
% parfor (ii = 1:length(subid), NSLOTS)
    %% Find all files with "ref#.mat" in the subfolder
    % Define entire filepath 
    fpath = fullfile(dpath, subid{ii}, subdir);
    % Create struct from directory contents
    list = dir(fpath);
    % Extract names
    names = {list.name};
    % Create regular expression
    if strcmp(subid{ii},'NC_7597')
        exp = 'Ref_BASIC\d*+.btf';
        % Filename to parse (this will be the same for each subject)
        fbase = 'Ref_BASIC';
    else
        exp = 'Ref_BASIC_4ds\d*+.btf';
        % Filename to parse (this will be the same for each subject)
        fbase = 'Ref_BASIC_4ds';
    end
    % Find strings matching regexp
    refnames = regexp(names, exp, 'match');
    % Remove empty elements (did not match)
    refnames = refnames(~cellfun('isempty',refnames));
    
    %% Combine the ref#.mat files into single matrix (c-scan)  
    % Load first ref stack for initializing matrix
    fname = strcat(fbase, '1.btf');
    filename = fullfile(fpath, fname);
    ref = TIFF2MAT(filename);    
    % Get dimensions of first stack
    [y,x,z] = size(ref);
    % Extrapolate total number of images in stack
    z = z .* length(refnames);
    % Initialize matrix for storing volume
    vol = zeros(y,x,z);
    % Add first reference to volume
    vol(:,:,1:size(ref,3)) = ref;
    % Indices to track last element index in z stack
    z0 = size(ref,3) + 1;
    % Iterate through ref#.mat files and add to vol
    for j=2:length(refnames)
        % Filepath to jth ref.mat file
        fname = strcat(fbase, num2str(j), '.btf');
        filename = fullfile(fpath, fname);
        % Load the ref stack
        tmp = load(filename);
        ref = tmp.Ref;
        % Add to volume matrix
        zf = z0+size(ref,3)-1;
        vol(:,:,z0:zf) = ref;
        z0 = zf + 1;
    end

    %% Save the output volume
    fout = fullfile(fpath, 'ref.mat');
    save_ref(fout, vol)
    % Save as tif
    segmat2tif(vol,fullfile(fpath, 'ref.tif'));
end

%% Function to save during parallelization
function save_ref(fout, vol)
% Save the stacked ref matrix
% INPUTS:
%   fout (string): the filepath for the output file to save
%   vol (double matrix): the stack of images to save

% Save the output
save(fout,'vol','-v7.3')

end







