clear 

addpath('/autofs/cluster/octdata2/users/Hui/tools/dg_utils/vol_recon_beta')
%% Enface 
nWorker = 12;
SaveTiff_Enface_script( ParameterFile, nWorker )
% SaveTiff_Enface_v2
% need to add orientation section
%% Create Macro.ijm file

WriteMacro_Spiral( ParameterFile );

%% Run ImageJ (FIJI) without ssh for speed. e.g. tigervnc

%% Run this to find median Fiji Tile Coordinates
%%%% (stored as Experiment_Fiji.X_Mean & .Y_Mean)
saveflag = 1;
[ Experiment_Fiji ] = ReadFijiTileCoordinates( ParameterFile, saveflag);

save(ParameterFile,'Experiment_Fiji','-append');
%% Mosaic 2D slices 

nWorker = 3;
Modality_2D = {'AIP','MIP','Retardance','Orientation'};%,'mus'};%Modality_2D = {'mus_minIP1', 'mus_AIP'}

for m = 1:4 %1:4
    modality = Modality_2D{m};       
    Mosaic2D_Telesto( ParameterFile, modality, nWorker );
end

%% StackNii
Modality_2D = {'AIP','MIP','Retardance','mus'};

for m = 1:3
    modality = Modality_2D{m};
    
    do_stacking(ParameterFile, modality)
end
%delete(findall(0,'type','figure','tag','TMWWaitbar'));

%% Mosaic 3D slices 

nWorker = 1;

modality = 'mus';       
Mosaic3D_Telesto( ParameterFile, modality, nWorker );

