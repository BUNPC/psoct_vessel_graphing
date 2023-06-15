function [ Experiment_Fiji ] = ReadFijiTileCoordinates( ParameterFile, saveflag )
% [ Experiment_Fiji] = ReadFijiRegisteredTileCoordinates( PathInfo, Acquisition, SliceInfo, FijiInfo, saveflag )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIAL SETTING 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(ParameterFile);

%%%% DISPLAY METHOD
Method      = ReadCoord.Method;

if      Method==1; fprintf(' -- Saving median coordinates -- \n');
elseif  Method==2; fprintf(' -- Saving reconstructed coordinates from median stepsize and offset -- \n');
end

%%%% INCLUDE SLICE
sliceid     = ReadCoord.ReadSliceID;
%%%% SAVE EXPERIMENT
fout        = ReadCoord.SaveExp;
fprintf(' - Save Experiment path = %s\n',fout);
%%%% LOAD TEMPLATE 
Experiment  = ReadCoord.ExperimentBasic;

%%%% SET COORDINATE & TILE BASENAME
if FijiInfo.isRegistered==1 
    fiji_file_ext = '.txt.registered';else;fiji_file_ext = '.txt';
end

pathname_coord = [FijiInfo.OutDir filesep FijiInfo.FileBase '%03i' fiji_file_ext];
pathname_tile  = [WriteMacro.TileDir filesep WriteMacro.Modality '_'];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Begin!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Experiment_Fiji = Experiment;

nslices = length(sliceid);
all_coords_x = zeros([size(Experiment_Fiji.X_Tot) nslices]);
all_coords_y = zeros([size(Experiment_Fiji.Y_Tot) nslices]);
MapIndex_Tot = zeros([size(Experiment_Fiji.MapIndex_Tot) nslices]);
for ssind=1:nslices
    currslice = sliceid(ssind);

    % Loop z tiles (should always have Z_tile==1 though)
    for zz = 1:Experiment_Fiji.Z_tile

%         % Get tile layout (tile numbers)
%         MapIndex = Experiment_Fiji.MapIndex_Tot(:,:,zz) + (currslice-1)*Experiment_Fiji.TilesPerSlice;
% 
%         % Get "Basic" tile coordinates (fixed overlap/uniform positions)
%         XX = squeeze(Experiment_Fiji.X_Tot(:,:,zz));
%         YY = squeeze(Experiment_Fiji.Y_Tot(:,:,zz));
% 
%         % Mtxs for Fiji optimal tile coordinates
%         X1 = XX; %zeros(size(MapIndex)); %XX;
%         Y1 = YY; %zeros(size(MapIndex)); %YY;
        
        
        
        MapIndex = Experiment_Fiji.MapIndex_Tot(:,:,zz)+(currslice-1)*Experiment_Fiji.TilesPerSlice;
        XX = squeeze(Experiment_Fiji.X_Tot(:,:,zz));
        YY = squeeze(Experiment_Fiji.Y_Tot(:,:,zz));
        X1 = zeros(size(MapIndex));
        Y1 = zeros(size(MapIndex));
    
    
    
    

        % Read Fiji stitching coordinates file
        %  (Can use 'registered' from compute_overlap, or simply
        %       TileCoordinates.txt from manual stitching in GUI)
        fid = fopen(sprintf(pathname_coord,currslice),'r');

        % Skip 1st 4 lines
        for ii=1:4
            kk = fgets(fid);
        end

        % Read coords
        ii=0;
        while ~feof(fid)
            ii = ii +1;
            kk = fgets(fid);
            img(:,ii)=sscanf(kk,[pathname_tile '%d.tiff; ; ( %f, %f)']);
        end
        fclose(fid);

        % Extract coords and save to matrices
        MapIndex2 = -1*ones(size(MapIndex));
        coord = (img');
        warningcounter = 0;
        for ii=1:size(coord,1)
            [row,col] = find(MapIndex==coord(ii,1));
            if isempty(row) && isempty(col) 
                warningcounter = warningcounter+1;
            end
            MapIndex2(row,col) = coord(ii,1);
            % round coords to nearest integer
            Y1(row,col) = round(coord(ii,2));
            X1(row,col) = round(coord(ii,3));
        end

        if warningcounter>0
            warning('Tile number mismatch detected between MapIndex and Fiji coord text file');
            disp(' -- Make sure ReadCoord.ReadSliceID is correct. ');
            disp(' -- Make sure tile nbs in MapIndex match those in text file.');
        end
       
    
        % Set min. coord value in X,Y = 1; and find nans;
        clear img coord;
        
        X        = X1;
        Y        = Y1;

        X           = X-min(min(X))+1;
        Y           = Y-min(min(Y))+1;
        indexNaN    = (MapIndex2==-1);
        
        
        X(indexNaN) = NaN;
        Y(indexNaN) = NaN;

        % Uniform the 'Origin' coordinate 
        X = X - X(round(size(X,1)/2),round(size(X,2)/2));
        Y = Y - Y(round(size(Y,1)/2),round(size(Y,2)/2));

        % Store coords in 3d matrix
        all_coords_x(:,:,ssind) = X;
        all_coords_y(:,:,ssind) = Y;
        MapIndex_Tot(:,:,ssind) = MapIndex2;
    end
end

Experiment_Fiji.X_Tot = all_coords_x;
Experiment_Fiji.Y_Tot = all_coords_y;
Experiment_Fiji.MapIndex_Tot = MapIndex_Tot;

% For case where Z_tiles>1
% Experiment_Fiji.X_Mean = squeeze(median(Experiment_Fiji.X_Tot(:,:,:),3));
% Experiment_Fiji.Y_Mean = squeeze(median(Experiment_Fiji.Y_Tot(:,:,:),3));
% Experiment_Fiji.X_std = squeeze(std(Experiment_Fiji.X_Tot(:,:,:),[],3));
% Experiment_Fiji.Y_std = squeeze(std(Experiment_Fiji.Y_Tot(:,:,:),[],3));

% Method 1: Set XY_Mean coords to the median of all slices (omit NaNs)
Method_1.X_Mean = round(squeeze(median(all_coords_x,3,'omitnan')));
Method_1.Y_Mean = round(squeeze(median(all_coords_y,3,'omitnan')));
Method_1.X_std = squeeze(std(all_coords_x,[],3,'omitnan'));
Method_1.Y_std = squeeze(std(all_coords_y,[],3,'omitnan'));
Experiment_Fiji.Method_1 = Method_1;

Experiment_Fiji.X_Mean = round(squeeze(median(all_coords_x,3,'omitnan')));
Experiment_Fiji.Y_Mean = round(squeeze(median(all_coords_y,3,'omitnan')));


% Method 2: Set XY_Mean coords from median offset and median step size
X_tile = Experiment_Fiji.X_tile;
Y_tile = Experiment_Fiji.Y_tile;

step_x          = median(   diff(all_coords_x, 1, 1),           'all','omitnan'); step_x = round(step_x);
matrix_x        = repmat(   [1: step_x : step_x * X_tile].',    [1 Y_tile]);
offset_X        = median(   diff(all_coords_x, 1, 2),           'all','omitnan'); offset_X = round(offset_X);
offset_map_X    = repmat(   [1:Y_tile]-1,                       [X_tile,1])         *offset_X;
Method_2.X_Mean = matrix_x + offset_map_X;
Method_2.step_x = step_x;
Method_2.offset_X = offset_X;

step_y          = median(   diff(all_coords_y, 1, 2),           'all','omitnan'); step_y = round(step_y);
matrix_y        = repmat(   [1: step_y : step_y * Y_tile],      [X_tile 1]);
offset_Y        = median(   diff(all_coords_y, 1, 1),           'all','omitnan'); offset_Y = round(offset_Y);
offset_map_Y    = repmat(   [1:X_tile].'-1,                     [1,Y_tile])         *offset_Y;
Method_2.Y_Mean = matrix_y + offset_map_Y;
Method_2.step_y = step_y;
Method_2.offset_Y = offset_Y;

Experiment_Fiji.Method_2 = Method_2;

switch Method
    case 1
        Experiment_Fiji.X_Mean = Experiment_Fiji.Method_1.X_Mean;
        Experiment_Fiji.Y_Mean = Experiment_Fiji.Method_1.Y_Mean;
        Experiment_Fiji.X_std = Experiment_Fiji.Method_1.X_std;
        Experiment_Fiji.Y_std = Experiment_Fiji.Method_1.Y_std;
    case 2
        Experiment_Fiji.X_Mean = Experiment_Fiji.Method_2.X_Mean;
        Experiment_Fiji.Y_Mean = Experiment_Fiji.Method_2.Y_Mean;
end


% Find tiles that were NaN for all slices - exclude in Mosaic)
Experiment_Fiji.Nantile_Prc = sum(isnan(Experiment_Fiji.X_Tot),3)/size(Experiment_Fiji.X_Tot,3);
% Experiment_Fiji.X_Mean(nantiles) = nan;
% Experiment_Fiji.Y_Mean(nantiles) = nan;
% 
% nantiles = Experiment_Fiji.Nantile_Prc>0.95;

nantiles = isnan(Experiment_Fiji.X_Mean);

Experiment_Fiji.MapIndex_Fiji = Experiment_Fiji.MapIndex_Tot_offset;
Experiment_Fiji.MapIndex_Fiji(nantiles) = -1;

fprintf('Y std = %.2f/%.2f (max/median [px]); \n', max(Experiment_Fiji.Y_std(:)),median(Experiment_Fiji.Y_std(:),'omitnan'));
fprintf('X std = %.2f/%.2f (max/median [px]); \n', max(Experiment_Fiji.X_std(:)),median(Experiment_Fiji.X_std(:),'omitnan'));


% Save Experiment_Fiji struct to .mat file (to load during Mosaic)
if saveflag==1
    save(fout,'Experiment_Fiji','FijiInfo','Experiment');
end








