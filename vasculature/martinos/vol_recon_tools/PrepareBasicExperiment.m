function [ ExperimentBasic ] = PrepareBasicExperiment( ParameterFile, TLSSLog, ExpBasic )

%
% [ ExperimentBasic ] = PrepareData2D_Basic_v2( ProcDir, file_TLSS_log, [ file_out ] )
%   Arguments:
%     Required-----
%       - ProcDir       = path to proc dir, containing OctScanInfoLog.txt
%       - file_TLSS_log = full path to TLSS log file for your scan
% 
%     (optional)
%       - file_out      = name of .mat file to save Experiment_Basic struct
%                           to. If not given, will not save file
% 
% 
%%% Experiment with 10x objective
%
% This script is used to create the files needed to run "Basic"
% mosaic/stitching of tiles for Octopus PSOCT acquisitions.
% 
% "Basic" refers to stitching with 50% overlap bewteen tiles, without using
% Fiji to "compute_overlap" using Fourier shift theorem. We have trouble
% using this Fiji approach for ex vivo human brain, likely due to lack of
% contrast/features in deep WM regions.
% 
% **** You can also specify any other value for percent overlap!!! ****
% 
% This also allows "prospective" data processing, so that you can visualize
% each stitched slice just a few minutes after it is acquired (almost real
% time!).
%
% This will create a structure ExperimentBasic, which contains fields
% X_mean and Y_mean. These are defining the pixel coordinates where the
% (top left) corner of each tile is placed within the mosaic/stitched
% slice, assuming 50% overlap b/w adjacent tiles. These are the parameters 
% you need to run Mosaic_Octopus_10xW_2D_Basic.m
% 
% You can run this immediately after starting the serial acquisition, then
% compile the Mosaic_Octopus_10xW_2D_Basic.m function, and run the shell
% script prospectiveDataProcessing.sh, which will run the stitching of each
% slice on the cluster, and send email updates once each slice is stitched.
% You can then monitor the status of the scan remotely, and see if there's
% a problem that requires you to go in and intervene/fix something.
%
%
% 
% TO RUN THIS SCRIPT, YOU NEED:
%  - name of TLSS log file from stage computer
%  - create an OctScanInfoLog.txt file inside ProcDir
% 
% 
% THIS SCRIPT CALLS THE FOLLOWING FUNCTIONS:
%  - ReadTLSSLogFile.m
% 
% 
% THIS SCRIPT DOES THE FOLLOWING:
%  - generates struct ExperimentBasic.mat
%  - if third argument is given, will save struct to the specified file
% 
%
%%%%%%%%%%%%%% EXAMPLE USAGE: %%%%%%%%%%%%%%
%   ProcDir = '/autofs/cluster/octdata2/users/Hui/Brainstem/processed/20211013_run2_processed';
%   file_TLSS_log = '/autofs/cluster/octdata2/users/Hui/Brainstem/processed/20211013_run2_processed/mosaic.txt';
%   file_out = [ProcDir filesep 'Experiment_Basic_run2.mat'];
%   [ ExperimentBasic ] = PrepareData2D_Basic_v2( ProcDir, file_TLSS_log, file_out );
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%    GET BASIC SCAN INFO FROM OctScanInfoLog.txt   %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(ParameterFile);
addpath('/autofs/cluster/octdata2/users/Hui/tools/rob_utils/OCTBasic/proc2d/tlss_utils');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%    SET SCAN-SPECIFIC PATHS+STUFF    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Load TLSS log file (from stage computer)
if ~exist(TLSSLog,'file')
    disp(['TLSS LOG FILE: ' TLSSLog]);
    error('logFile:doesNotExist',...
        'Log file does not exist!');
end
fid = importdata(TLSSLog,'\n'); 

%Can just use only first slice, no need to do all slices
sliceid = 1;

for ss = 1:length(sliceid)
    disp(['slice ' num2str(sliceid(ss))]);
    
    %%%% To read in experiment info (from stage computer)
    disp(' Reading TLSS log file...');
    ExperimentBasic = ReadTLSSLogFile(fid, sliceid(ss));
    
    %%% If telesto, flip mapindex (tiles start in top right, not top left)
    if strcmpi(Scan.System,'Telesto')
        ExperimentBasic.MapIndex_Tot = fliplr(ExperimentBasic.MapIndex_Tot);
    end
    
    %%%% To set tile FOV (microns)
    ExperimentBasic.FOV = Scan.FoV;  %in micron
    %%%% To set size of fov in pixels in x,y
    ExperimentBasic.NbPix = Scan.NbPixels; 
    %%%% To set first tile number (in case offset in file numbering)
    ExperimentBasic.First_Tile = Scan.First_Tile;
    ExperimentBasic.MapIndex_Tot_offset = ExperimentBasic.MapIndex_Tot + ExperimentBasic.First_Tile -1;
    %%%% To set % overlap between tiles in x,y
    ExperimentBasic.PercentOverlap = round((1-Scan.StepSize/Scan.NbPixels)*100);
    %%%% To calculate pixels size (microns)
    ExperimentBasic.PixSize = ExperimentBasic.FOV / ExperimentBasic.NbPix;
    
    
    %%%% To calculate pixels per step in x,y (px between adjacent tiles)
    % way 1 - using [nbpix in x,y] and [percent overlap in x,y]
%     ExperimentBasic.X_step_pix = ExperimentBasic.NbPix * (ExperimentBasic.PercentOverlap * 1e-2);
%     ExperimentBasic.Y_step_pix = ExperimentBasic.NbPix * (ExperimentBasic.PercentOverlap * 1e-2);
%     ExperimentBasic.Z_step_pix = 1;
    % way 2 - using [step size microns in x,y] and [pixel size in x,y]
    ExperimentBasic.X_step_pix = ExperimentBasic.X_step / ExperimentBasic.PixSize;
    ExperimentBasic.Y_step_pix = ExperimentBasic.Y_step / ExperimentBasic.PixSize;


    %%%% To convert X,Y_Tot from micron to pixels 
    ExperimentBasic.X_Tot_micron = ExperimentBasic.X_Tot;
    ExperimentBasic.Y_Tot_micron = ExperimentBasic.Y_Tot;
    
    ExperimentBasic.X_Tot = round(ExperimentBasic.X_Tot_micron / ExperimentBasic.PixSize);
    ExperimentBasic.Y_Tot = round(ExperimentBasic.Y_Tot_micron / ExperimentBasic.PixSize);
    
    ExperimentBasic.X_Tot = ExperimentBasic.X_Tot - min(ExperimentBasic.X_Tot(:)) + 1;
    ExperimentBasic.Y_Tot = ExperimentBasic.Y_Tot - min(ExperimentBasic.Y_Tot(:)) + 1;
    indexNaN = (ExperimentBasic.MapIndex_Tot==-1);
    ExperimentBasic.X_Tot(indexNaN)=round(ExperimentBasic.X_Tot(indexNaN) - min(ExperimentBasic.X_Tot(:)) + 1);
    ExperimentBasic.Y_Tot(indexNaN)=round(ExperimentBasic.Y_Tot(indexNaN) - min(ExperimentBasic.Y_Tot(:)) + 1);
    
    %%%% Save as *_Mean field, as with previous notation
    ExperimentBasic.X_Mean = squeeze(median(ExperimentBasic.X_Tot,3));
    ExperimentBasic.Y_Mean = squeeze(median(ExperimentBasic.Y_Tot,3));
    ExperimentBasic.X_std = squeeze(std(ExperimentBasic.X_Tot,[],3));
    ExperimentBasic.Y_std = squeeze(std(ExperimentBasic.Y_Tot,[],3));

    %%%% To run Mosaic with 50% overlap, call ExperimentBasic.*_Mean to get
    %%%%   the x and y pixel positions of each tile in the slice.
end

if ~isempty(ExpBasic)
    %%%% To save ExperimentBasic struct to .mat file
    fprintf(' - Writing ExperimentBasic struct to .mat file:\n %s\n',ExpBasic);
    
    save(ExpBasic,'ExperimentBasic');
end

end
