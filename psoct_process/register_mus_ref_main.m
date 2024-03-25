%% Create a mask for the white matter
% Use the scattering map to segment just the white matter
%{
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

%% Initialize parallel pool
% Retrieve the number of available cores
n_cores = str2num(getenv('NSLOTS'));
% Set the maximum number of threads equal to the number of cores
maxNumCompThreads(n_cores);

% Check whether parallel pool exists
% poolobj = gcp('nocreate');
% if isempty(poolobj)
%     % No parallel pool exists. Initialize a parallel pool.
%     pc = parcluster('local');
%     % Setup directory for logging
%     pc.JobStorageLocation = getenv('TMPDIR');
%     % Start running the parallel pool
%     parpool(pc, n_cores);
% end

%% Initialize subject ID lists
%%% All subjects to analyze
subid = {'AD_10382', 'AD_20832', 'AD_20969',...
         'AD_21354', 'AD_21424',...
         'CTE_6489', 'CTE_6912',...
         'CTE_7019', 'CTE_7126',...
         'NC_6839',  'NC_6974', 'NC_8653',...
         'NC_21499', 'NC_301181'};

%%% Number of scattering maps (mus_#.tif) per volume
% Note, the scattering maps for CTE_8572 must be generated
nmus = [17,15,22,15,20,...
        22,10,18,18,...
        16,20,15,15,24];
% Number of frames in z-dimension per physical slice
nstack = 11;
% Number of frames in z-dimension per volume. This is dependent on the
% number of mus maps that were processed for each tissue volume.
zmax = nmus .* nstack;

% Create dictionary to store last image in stack
regstruct = struct;
for ii = 1:length(subid)
    sub = char(subid(ii));
    regstruct.(sub).zmax = zmax(ii);
end

% Debugging variable. True to view figures.
viz = true;

%% Initialize directories and filenames

%%% Directories and file name
% Upper level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T';    
% Subfolder containing ref files
subdir = '/dist_corrected/volume/ref';
% Scattering subdirectory
scatdir = 'fitting_4x/';
% non-normalized PSOCT volume
vol_name = 'ref';

%% Iterate over subjects.
for ii = 1:length(subid)
    %%% Debugging information
    fprintf('\n---------Starting Subject %s---------\n',subid{ii})

    %%% Import non-normalized volume
    fprintf('importing non-normalized volume\n')
    fpath = fullfile(dpath, subid{ii}, subdir, strcat(vol_name,'.mat'));
    vol = load(fpath); 
    vol = vol.vol;

    %%% Retrieve last slice # and N slices per physical stack
    sid = char(subid{ii});
    % Retrieve the last slice #
    zmax = regstruct.(sid).zmax;
    
    %% Iterate over frames of OCT intensity volume
    % Initialize array of z-stack frames that mark a new physical slice.
    % This indicates that the next frame in the mus stack should be taken.
    z = 1:nstack:zmax;
    % Counter for the number of mus frames
    mus_cnt = 1;
    % Init. matrices for registration stats and registered scattering map
    dfield_mat = zeros(size(vol,1),size(vol,2),size(vol,3).*2);
    mus_mat = zeros(size(vol));
    
    % Iterate over slices
    for j = 1:zmax
        %%% Debugging information
        fprintf('Starting Frame # %i \n', j)

        %% Import scattering Map
        % Import the scattering map for each physical slice
        if ismember(j,z)
            fprintf('importing mask\n')
            % Create string of scattering map
            mus_str = strcat('mus',num2str(mus_cnt),'.tif');
            fpath = fullfile(dpath, sid, scatdir, mus_str);
            % Rescale the scattering map to uint8 for registration
            mus = uint8(rescale(TIFF2MAT(fpath),0,255));
            % Iterate the mus slice count
            mus_cnt = mus_cnt + 1;
        end
        %% Register scattering map to first slice
        % Overlay pair prior to registration
        slice = vol(:,:,j);
        slice = uint8(rescale(slice,0,255));
        % Perform registration
        mus_reg = register_mus_ref_slice(mus, slice);
        dfield = mus_reg.DisplacementField;
        mus_reg = mus_reg.RegisteredImage;
        
        %%% Display before/after registration
        if viz
            % Full screen figure
            figure('units','normalized','outerposition',[0 0 1 1]);
            % Before registration
            subplot(1,2,1); imshowpair(slice, mus); title('Unregistered')
            % Show the pair of images
            subplot(1,2,2); imshowpair(slice, mus_reg); title('Registered')
            %%% Save figure
            % Check if subfolder exists
            regdir = fullfile(dpath, sid, scatdir, 'registered');
            if ~exist(regdir,'dir')
                mkdir(regdir)
            end
            % Create file path to output registered images
            fname = strcat(mus_str(1:end-4),'_slice',num2str(j),'_reg.png');
            reg_plot = fullfile(dpath, sid, scatdir, 'registered',fname);
            saveas(gcf,reg_plot);
            close all;
        end

        %%% Store the registration and registered image
        mus_mat(:,:,j) = mus_reg;
        dfield_mat(:,:,j:j+1) = dfield;
    end
    %%% Store the masks and masked tissue in structure
    regstruct.(sid).dfield = dfield_mat;
    regstruct.(sid).mus_reg = mus_mat;
    % Save the struct
    fout = fullfile(dpath, sid, scatdir, 'registered','regstruct.mat');
    save(fout,'regstruct');

    %%% Debugging Info
    fprintf('\n---------Finished Subject %s---------\n',subid{ii})
end

