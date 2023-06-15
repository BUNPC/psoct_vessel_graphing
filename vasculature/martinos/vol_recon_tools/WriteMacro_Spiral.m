function WriteMacro_Spiral( ParameterFile ) 

% This script delet previous Macro
% This script consider Agar info from script SaveTiff_Enface_script
% This script write out ImageJ reading sequence in a spiral-out order from
% a set origin

% fdir = '/autofs/cluster/octdata2/users/Hui/SLF_I55/slab4/ProcessI55_slab4_x0z1_20210828/Stitching';
% fmacro = 'Fiji_Stitching_Macro_clipped.ijm';
% SliceInfo.sliceid = 90;
% SliceInfo.modality = 'AIP';
% SliceInfo.tiledir = '/autofs/cluster/octdata2/users/Hui/SLF_I55/slab4/ProcessI55_slab4_x0z1_20210828/Enface_Tiff_clipped';
% load('/autofs/cluster/octdata2/users/Hui/SLF_I55/slab4/ProcessI55_slab4_x0z1_20210828/StitchingBasic/ExperimentBasic.mat');
% Experiment = ExperimentBasic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIAL SETTING 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(ParameterFile);

%%%% CHECKING MACRO
fmacro      = WriteMacro.SaveMacro;
if exist(fmacro,'file'); delete(fmacro);
fprintf(' --deleting existing Macro file-- \n%s\n',fmacro);end
%%%% SETTING PARAMETERS
origin      = WriteMacro.SpiralOrigin;
modality    = WriteMacro.Modality;
sliceid     = WriteMacro.SliceID;
%%%% SETTING DIRECTORIES 
tiledir     = WriteMacro.TileDir;
AgarInfoDir = WriteMacro.AgarInfoDir;
fdir        = WriteMacro.FijiOutDir; if ~exist(fdir,'dir'); mkdir(fdir); end

fprintf(' - Tile directory = \n%s\n',tiledir);
fprintf(' - Fiji Input/Output coordinate directory = \n%s\n',fdir);
fprintf(' - Agar directory = \n%s\n',AgarInfoDir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING MOSAIC PARAMETERS 
fprintf(' - Loading Experiment file...\n %s\n', Parameters.ExpBasic);
S = whos('-file',Parameters.ExpBasic);
load(Parameters.ExpBasic);
Experiment  = eval(S(find(contains({S(:).name},'Experiment'))).name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Begin!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid_Macro = fopen(fmacro, 'a'); %'a');

for ssind=1:length(sliceid)
    currslice = sliceid(ssind);
    slicepad = sprintf('%03i',currslice);
    
    for zz = 1:Experiment.Z_tile
        X = (round(squeeze(Experiment.X_Tot(:,:,zz))));
        Y = (round(squeeze(Experiment.Y_Tot(:,:,zz))));
        MapIndex = Experiment.MapIndex_Tot(:,:,zz) + (currslice-1)*Experiment.TilesPerSlice;

        % Create the file for the individual stitching of each depth
        fmosaic = [fdir filesep 'Mosaic_depth' sprintf('%03i',zz) '_slice' slicepad '.txt'];
        fid = fopen(fmosaic, 'w');

        l1='# Define the number of dimensions we are working on\n';
        fprintf(fid,l1);
        l2='dim = 2\n \n';
        fprintf(fid,l2);
        l3='# Define the image coordinates\n';
        fprintf(fid,l3);  

        
        %spiral out
        mapsz = size(MapIndex);
        
        grid_start = 50 - origin + 1;
        grid_end = grid_start + mapsz -1;
        B = ReadSpiral(99,99);
        B=fliplr(B);
%        if flipflag ==1; B=fliplr(B);end
        B = B(grid_start(1):grid_end(1),...
            grid_start(2):grid_end(2));

        
        [~,id] = sort(B(:),'descend');
        for ii = 1:length(id)
            tilepad = sprintf('%03i', MapIndex(id(ii)));
            load([AgarInfoDir filesep tilepad '.mat'], 'tileinfo');
            if tileinfo(2)==0 % if tileinfo(2) is 0 which is agar, then skip
                continue
            end
            if tileinfo(1)~=MapIndex(id(ii))
                fprintf('current tile %s and mapindex tile %i doesnt match', tilepad, tileinfo(1))
            end
            
            lig=sprintf('%s/%s_%03i.tiff; ; (%03i,%03i)\n', tiledir,modality,MapIndex(id(ii)),Y(id(ii)),X(id(ii)));
%             lig=[pathname_tile tilepad '.tiff; ; (' sprintf('%03i',Y(id(ii))) ',' sprintf('%03i',X(id(ii))) ')\n'];
            fprintf(fid,lig);
        end
        
        fclose(fid);

        % Create the command line for the actual stitching
        m1=['run("Stitch Collection of Images", "layout=' fmosaic ' compute_overlap channels_for_registration=[Red, Green and Blue] rgb_order=rgb fusion_method=[Linear Blending] fusion=1.50 regression=0.30 max/avg=2.50 absolute=3.50");\n'];
        %m2=['saveAs("Tiff", "' fdir filesep fspiral filesep SliceInfo.modality '_Slice' slicepad '.tif");\n'];
        m2=['saveAs("Tiff", "' fdir filesep modality '_Slice' slicepad '.tif");\n'];
        m3 = 'close();\n';
        fprintf(fid_Macro,m1);
        fprintf(fid_Macro,m2);
        fprintf(fid_Macro,m3);

    end
end

fclose(fid_Macro);
