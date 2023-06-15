clear
addpath('/autofs/cluster/octdata2/users/Hui/tools/dg_utils/vol_recon_beta');
logflag = false;
ParameterFile  = ['/autofs/cluster/octdata2/users/Chao/caa/caa_6/occipital/process_20220209_run3/Parameters.mat'];
%P = whos('-file',ParameterFile);
if exist(ParameterFile,'file'); return; end
%% Aquisition Parameters 
Scan.SampleName          = 'CAA_6 Occipital Run3';
Scan.Date                = '20220106';
Scan.OCTOperator         = 'Bill';
Scan.PostProcessOperator = 'Dylan';
Scan.ImagingNotes        = 'Saved crop volume, saved surface, saved complex...';
Scan.FilePrefix          = 'test_processed_';
Scan.Objective           = '40x';
Scan.System              = 'Telesto'; IsTelesto = strcmpi(Scan.System,'Telesto');

Scan.RawDataDir          = '/autofs/';
Scan.ProcessDir          = '';
Scan.TLSS_log            = 'run3.txt'; % needed?
Scan.First_Tile          = 1;

Scan.Thickness           = 50;   % um  50um serial scanning; 100um cutting
Scan.FoV                 = 3600;  % um
Scan.NbPixels            = 360;   % pix
Scan.StepSize            = 320;   % pix
Scan.CropDepth           = 100;   % pix  size(niftiread([Scan.RawDataDir '/test_processed_001_cropped.nii']));

in_planeRes              = 0.010;  % mm % Scan.FoV/Scan.NbPixels/1000;
thru_planeRes            = 0.0025;% mm
Scan.Resolution          = [in_planeRes in_planeRes thru_planeRes];

save(ParameterFile,'Scan');
%% Folder & File names
ProcDir     = Scan.ProcessDir;
RawDataDir  = Scan.RawDataDir;

Enface_Tiff    = [ProcDir filesep 'Enface_Tiff'];
Enface_MAT     = [ProcDir filesep 'Enface_MAT'];
StitchingFiji  = [ProcDir filesep 'StitchingFiji'];
StitchingBasic = [ProcDir filesep 'StitchingBasic'];
StackNii       = [ProcDir filesep 'StackNii'];
MapIndex       = [ProcDir filesep 'MapIndex'];

OctScanInfoLog = [ProcDir filesep 'OctScanInfoLog.txt'];
TLSSLog        = [ProcDir filesep 'run3.txt'];
ExpBasic       = [ProcDir filesep 'ExperimentBasic.mat'];
ExpFiji        = [StitchingFiji filesep 'Experiment_Fiji.mat'];
% Macro_ijm is defined later.

%% ExperimentBasic
if ~exist(ExpBasic,'file')
    [ ExperimentBasic ] = PrepareBasicExperiment( ParameterFile, TLSSLog, ExpBasic );
    save(ParameterFile,'ExperimentBasic','-append');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
load(ExpBasic);
Parameters.SliceID = 1:5;
Parameters.TileID  = 1:ExperimentBasic.TilesPerSlice * Parameters.SliceID(end);

Parameters.OrientationSign = 1;
Parameters.OrienOffset     = 0;

Parameters.XPixClip = 18;
Parameters.YPixClip = 0;

Parameters.ExpFiji = ExpFiji;
Parameters.ExpBasic= ExpBasic;
%Mosaic2D.InFileType = 'nifti';
save(ParameterFile,'Parameters','-append');

%%
if ~exist('Parameters','var');load(ParameterFile);end

% Parameters.AipGrayRange = 
% Parameters.MipGrayRange = 
% Parameters.RetGrayRange = 

GrayRangeCheck = ...
   [~isfield(Parameters,'AipGrayRange'),...
    ~isfield(Parameters,'MipGrayRange'),...
    ~isfield(Parameters,'RetGrayRange')];%, ...isfield(Parameters,'musGrayRange'); 

if any(GrayRangeCheck)
    modalid = find(GrayRangeCheck);

    Mosaic2D.sliceidx   = [1;1;1]; %slice2 in; slice2 out; run1
    Mosaic2D.imresize   = 0; 

    Mosaic2D.InFileType = 'nifti'; % 'mat'
    Mosaic2D.indir      = {RawDataDir}; % Enface_MAT
    Mosaic2D.outdir     = StitchingBasic;
    Mosaic2D.Exp        = ExpBasic;
    save(ParameterFile,'Mosaic2D','-append');

    % Find GrayRange - Preparation for Fiji Registration
    adju = [5,1;12,1;3,-5];
    modals = {'AIP','MIP','Retardance','mus'};
    modalsstr = {'AipGrayRange', 'MipGrayRange', 'RetGrayRange','musGrayRange'};
    for i = modalid
        modality = modals{i};
        Mosaic2D_Telesto( ParameterFile, modality, 1 );
        load([StitchingBasic filesep modality '_slice001.mat'],'MosaicFinal');
        GrayRange(1) = floor(prctile(MosaicFinal(MosaicFinal~=0),0.1))+adju(i,1); 
        GrayRange(2) = ceil(prctile(MosaicFinal(MosaicFinal~=0),99.9))+adju(i,2);
        Parameters.(modalsstr{i}) = GrayRange;
        %figure; imagesc(MosaicFinal); axis tight off equal; caxis([Parameters.(modalsstr{modalid})]); colormap(gray); title([modality ' ' num2str(Parameters.(modalsstr{modalid})(1)) ' ' num2str(Parameters.(modalsstr{modalid})(2))]);
    end
    save(ParameterFile,'Parameters','-append');
end
%% 
Parameters.Agar.threshold_pct = 0.04;
Parameters.Agar.threshold_std = 1.5; % will this change in TDE?
Parameters.Agar.gm_aip = diff(Parameters.AipGrayRange)/4+Parameters.AipGrayRange(1);%(find thresh at mid-range intensity of graymatter which is roughly at lower half of AipGrayRange)
Parameters.Agar.dir    = MapIndex;
save(ParameterFile,'Parameters','-append');

% use std on mat2gray(I,aipgrayscale) may yeild a more distinct result, std
% thresh may have variation then

%[ agar_status ] = view_tile_aip( ParameterFile, [7 5] );
% SaveTiff_Enface_script( ParameterFile )

Enface.indir  = RawDataDir;
Enface.outdir = Enface_Tiff;
%Enface.save = {'AIP','MIP','Retardance','Orientation','mus'};
Enface.save = {'AIP','MIP','Retardance'};%,'Orientation','mus'};
save(ParameterFile,'Enface','-append');
%% WriteMacro & Read Registered Coordinates
WriteMacro.Modality    = 'AIP'; % 'AIP' 'mus'
WriteMacro.SpiralOrigin= round([ExperimentBasic.X_tile,ExperimentBasic.Y_tile]/2);

WriteMacro.SliceID     = Parameters.SliceID;
WriteMacro.TileDir     = Enface_Tiff;
WriteMacro.FijiOutDir  = StitchingFiji;
WriteMacro.SaveMacro   = [StitchingFiji filesep 'Macro_' WriteMacro.Modality '.ijm'];
WriteMacro.AgarInfoDir = MapIndex;
save(ParameterFile,'WriteMacro','-append'); % FijiStep

ReadCoord.Method       = 1;     %select a method. 1=median; 2=step&offset
ReadCoord.ReadSliceID  = Parameters.SliceID;
ReadCoord.CoordDir     = StitchingFiji;
ReadCoord.SaveExp      = ExpFiji;
load(ExpBasic, 'ExperimentBasic');
ReadCoord.ExperimentBasic=ExperimentBasic;

FijiInfo.isRegistered  = 1;
FijiInfo.OutDir        = StitchingFiji;
FijiInfo.FileBase      = 'Mosaic_depth001_slice';
save(ParameterFile,'FijiInfo', 'ReadCoord', '-append');

% EXAMPLE: notes: Fiji tif output jumping around
% EXAMPLE: Y std = 57.81/9.55 (max/median); 
% EXAMPLE: X std = 14.53/9.34 (max/median);
%% 
tmp.in  = [Parameters.SliceID];                                                      %sliceid from multiple runs
tmp.out = 1:length(tmp.in);             %
tmp.run = [ones(size(Parameters.SliceID))*1];
sliceidx = [tmp.in;tmp.out;tmp.run];

%% Matlab Mosaic Parameter
Mosaic2D.sliceidx   = sliceidx;
Mosaic2D.imresize   = 0; 

Mosaic2D.InFileType = 'nifti'; % 'mat'
Mosaic2D.indir      = {RawDataDir}; % Enface_MAT
Mosaic2D.outdir     = StitchingFiji;
Mosaic2D.Exp        = ExpFiji;

if logflag == 1; RecordHistory(ParameterFile);end
save(ParameterFile,'Mosaic2D','-append');

%% StackNii
Stack.indir   = StitchingFiji;
Stack.outdir  = StackNii;
Stack.sliceid = sliceidx(2,:);
Stack.Resolution = [Scan.Resolution(1) Scan.Resolution(2) Scan.Thickness/1000];
save(ParameterFile,'Stack','-append');

%% Mus
% SaveMus.depth = 100;
% SaveMus.z_mean= 5:80;
% SaveMus.save %2d tile 3d tile%
% SaveMus.Method %Vermeer % Optical_fitting_MGH 

%% FijiParameters
% channels_for_registration=[Red, Green and Blue] 
% rgb_order=rgb 
% fusion_method=[Linear Blending] 
% fusion=1.50 
% regression=0.30 
% max/avg=2.50 
% absolute=3.50

%% 

% rawdatadir
%
% ['test_processed_' sprintf('%03i',i) '_' lower(modality) '.nii']
% sprintf('test_processed_%03i_%s.nii',i,lower(modality));
% test_processed_001_aip.nii
% test_processed_001_mip.nii
% test_processed_001_orientation.nii
% test_processed_001_retardance.nii

% test_processed_001_cropped.nii
% test_processed_001_surface_finding.nii


% Enface_Tiff
%
% AIP_001.tiff
% MIP_001.tiff
% Retardance_001.tiff
% Orientation1_001.tiff
% Orientation2_001.tiff


% StitchingFiji
% 
% AIP_slice001.mat
% MIP_slice001.mat
% Retardance_slice001.mat
% Orientation_slice001.mat
% mus_slice001.mat
% 
% AIP_slice001.tiff
% AIP_Slice001.tif                          <- fiji output
% MIP_slice001.tiff
% Retardance_slice001.tiff
% Orientation1_slice001.tiff
% Orientation2_slice001.tiff
% 
% resized-AIP_slice001.jpg
% resized-MIP_slice001.jpg
% resized-Retardance_slice001.jpg
% resized,nearest-Orientation1_slice001.tiff
% 
% Mosaic_depth001_slice001.txt              <- WriteMacro output
% Mosaic_depth001_slice001.txt.registered   <- fiji output


% StackNii
% 
% Stacked_AIP.nii.gz 
% Stacked_MIP.nii.gz  
% Stacked_Retardance.nii.gz 
% Stacked_mus.nii.gz


