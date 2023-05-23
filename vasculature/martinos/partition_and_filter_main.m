%{
The purpose of this script is to divde the volume into smaller volumes.
This is to ensure that a computer has sufficient memory for applying
the Frangi filter.

TODO:
- this is currently breaking on line 166 of MRIread.m because the header is
      empty.
- programatically determine mat2cell variables
- programatically change hte skirt size.
	- skirt size = (# in conv3D(I,#) in Frangi filter / 2)
	- line "I = conv3D(I,40)" in "frangi_ves4_..m"
- Create function for running BASH:
	- change both input and output names

%% create_mri_mosaic.py commandline 
%  create_mri_mosaic.py is Jackson's script to divide large nii image into 
%  small blocks using freesurfer functions. The small blocks retained its
%  *original positions* when loaded in freeview. 

%% create_mri_mosaic.py -> I_seg (frangi filter output)
% First line is the function
% Second line is the input data.
% [1100 1000 1200] is the dimension of each block. This needs to be adjusted manually to match the dataset.
% [3 3 1] sets the number of subdivisions of the volume. This is to reduce the filesize.
%	Here, a 50GB file will be divided into 9 volumes.
/autofs/space/tiamat_001/users/jackson/projects/scripts/Python/create_mri_mosaic.py \ 
/autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/process_run1/mus_vessel/mus_mean_10um-iso.slice40px.I_seg.nii 1100 1000 1200 --m_dims 3 3 1   -s --outdir /autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/process_run1/mus_vessel/nii_split --tile_name seg

%% create_mri_mosaic.py -> I (PSOCT image of mus)
% The intern will use the PSCOT image generated here as the background.
% Then, they will overlay the segmentation.
% 
% First line is the function
% Second line is the input data.
/autofs/space/tiamat_001/users/jackson/projects/scripts/Python/create_mri_mosaic.py \
/autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/process_run1/mus_vessel/mus_mean_10um-iso.slice40px.t3_9.sigma1.nii 1100 1000 1200 --m_dims 3 3 1   -s --outdir /autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/process_run1/mus_vessel/nii_split --tile_name I

%}

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

%% Load output of stitching
fname = 'mus_mean_10um-iso.slice40px.nii';
fout = 'mus_mean_10um-iso.slice40px.I_seg.nii';

if ispc
    dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\CAA25_Frontal\Process_caa25_frontal\mus_vessel\';
else
    dpath = '/projectnb/npbssmic/ns/CAA25_Frontal/Process_caa25_frontal/mus_vessel/';
end

% Create path to save partitions in same folder as original data
partition_path = strcat(dpath, 'Ma_partition_cell_padded.mat');

%{
mri = MRIread(strcat(dpath, fname));
fprintf('mri loaded.\n')
mri.vol = single(mri.vol);

%% Divide volume into smaller sections for memory management.
% Programatically select these values to generate 18 volumes.
% (i.e.  divide entire volume into (3x3x2))

cell_I = mat2cell(mri.vol,[1000,1000,879],[1000,1000,613],[300, 390]);
clear mri
%% Initialize padded partitions
cell_I_padded = cell(size(cell_I)+2);
cell_I_padded = cellfun(@single,cell_I_padded,'UniformOutput',false);
cell_I_padded2 = cell_I_padded;
cell_I_m_padded = cell_I_padded;
cell_I_padded(2:end-1,2:end-1,2:end-1) = cell_I; % save partitioned images in the center

%% Set skirt size for padding subvolumes.
% The skirt size = (# in conv3D(I,#) in Frangi filter / 2) (I = conv3D(I,40))
% the padding size [#] can be set with different value.
% the skirt of the image which contains artifact will be trimmed after 
% the image processing and the frangi filter.
str = ["end-20+1:end",":","1:20"];
[x,y,z]=meshgrid(1:3);
str_idx = [y(:),x(:),z(:)];

% set linear index for a 3*3*3 moving cube window
lnridx_neighoor = sub2ind(size(cell_I_padded),y(:),x(:),z(:));
lnridx_neighoor = lnridx_neighoor - lnridx_neighoor(14); % off set index to be relative to the center block (14th)

%% add skirt to image
lnridx = find(cellfun(@isempty,cell_I_padded)==0);

for i=1:length(lnridx)
    % extract a 3*3 cube from cell_I_padded
    % cell_cube stores image detail; cell_cube_mask is binary mask of skirt
    cell_cube = cell_I_padded(lnridx(i)+lnridx_neighoor);
    cell_cube_mask = cellfun(@true_like,cell_cube,'UniformOutput',false);
    cell_cube_mask{14}(:) = false;

    % trim neighoor image block into skirt in cell_cube
    for ii = 1:length(str_idx)
        if isempty(cell_cube{ii})
            continue
        end
        str_ = sprintf('(%s,%s,%s)',str(str_idx(ii,:)));
        eval(['cell_cube{ii} = cell_cube{ii}' str_ ';'])
        eval(['cell_cube_mask{ii} = cell_cube_mask{ii}' str_ ';'])
    end
    
    cell_cube = reshape(cell_cube,3,3,3);
    % cut away empty cell (none-skirt and none-mainbody block)
    cell_cube = cut_empty_cell(cell_cube);
    % sew together the skirt with mainbody
    cell_I_padded2{lnridx(i)} = cell2mat(cell_cube);

    cell_cube_mask = reshape(cell_cube_mask,3,3,3);
    cell_cube_mask = cut_empty_cell(cell_cube_mask);
    cell_I_m_padded{lnridx(i)} = cell2mat(cell_cube_mask);
end

cell_I2 = cell_I_padded2(2:end-1,2:end-1,2:end-1);                          % padded image in each partition slot
cell_I2_m = cell_I_m_padded(2:end-1,2:end-1,2:end-1);                       % skirt mask of padded image in each partition slot
cell_I2_skirt_lnridx = cellfun(@find,cell_I2_m,'UniformOutput',false);      % linear index of skirt (save storage space)
cell_I2_skirt_lnridx = cellfun(@uint32,cell_I2_skirt_lnridx,'UniformOutput',false); % (save storage space)
cell_I2_main = cellfun(@size,cell_I,'UniformOutput',false);                 % save orignal partitioned block size for reconstruction
cell_I2_main_skirt = cellfun(@size,cell_I2,'UniformOutput',false);          % save padded partitioned block size
clear cell_I_padded* cell_I2_m


%% save partition
% Initialize output variable
placeholder = 1;
save(partition_path, 'placeholder','-v7.3')

for s = 1:length(cell_I2(:))
    var_name = sprintf('I%i',s);
    eval([var_name ' = cell_I2{s};'])
    save(partition_path, var_name,'-append')
    eval(['clear ' var_name])
end
reshape_sz = size(cell_I2);

% Save variables in to partition output structure
save(partition_path,'reshape_sz','-append')
save(partition_path,'cell_I2_skirt_lnridx','-append')
save(partition_path,'cell_I2_main','-append')
save(partition_path,'cell_I2_main_skirt','-append')

%}
%% Image Processing and Frangi Filter
frangi_partitions(dpath, partition_path);

%% I_seg reconstruction
load(partition_path,'reshape_sz','cell_I2_main','cell_I2_main_skirt','cell_I2_skirt_lnridx');

cell_seg = cell(reshape_sz);

% Remove padding from each tile
for s = 1:length(cell_I2_main(:))
    
    var_in_name = sprintf('I_seg%i',s);
    
    load(strcat(dpath,'ves4_padded.mat'), var_in_name)
    eval(['I_seg = ' var_in_name ';'])
    eval(['clear ' var_in_name])

    I_seg(cell_I2_skirt_lnridx{s})=[];
    I_seg = reshape(I_seg, cell_I2_main{s});
    cell_seg{s} = I_seg;
end
I_seg_all = cell2mat(cell_seg);

% Save to same directory as original data
save_mri_s(I_seg_all, strcat(dpath, fout), [0.01,0.01,0.01], 'uchar')

%% Smooth volume, convert to gray, normalize, Frangi filter
function frangi_partitions(dpath, ppath)
% INPUT:
%       ppath (string): path to the partitions

% Clear memory and reload partitions
load(ppath,'cell_I2_main');

% Define window size
window = [3 19];

% Save output variable

save(strcat(dpath, 'ves4_padded.mat'),'window','-v7.3')

for s = 1:length(cell_I2_main(:))
    tic
    var_in_name = sprintf('I%i',s);
    load(ppath,var_in_name)
    eval(['I = ' var_in_name ';'])
    eval(['clear ' var_in_name])
    
    I = smooth3(I,'gaussian',9,5);  % smooth sigma = 5
    I = mat2gray(I,[3 19]);         % window I [3 19]
    I = conv3D(I,40);               % normalize, kernel 40px
    I = mat2gray(I,[0 1]);          % window I [0 1]
    I = double(I);                  % make I double before frangi
    I_seg = ves(I);                 % frangi sigma [1 3 5]; I_VE thresh = 0.15; bwconncomp thresh = 100;                           
    toc
    tic
    var_out_name = sprintf('I_seg%i',s);
    eval([var_out_name ' = I_seg;'])
    clear I_seg
    save(strcat(dpath,'ves4_padded.mat'),var_out_name,'-append')
    eval(['clear ' var_in_name])
    toc
end

end

%% Frangi filter wrapper
function I_seg = ves(I)
    
    sigmas = [ 1,   3,    5]; %% 100~200 um 20px big vessel are manually segmented
    thres  = 0.15; 
    cthres = 100;

    sz = size(sigmas);
    I_seg = zeros(size(I),'logical');
    for i = 1:sz(2)
        tic
        I_VE = vesSegment(1-I,[sigmas(i)],1);   % I_VE is logical 
        I_VE = I_VE > thres;
        I_VE = clean_ves(I_VE,cthres);
        toc        
        I_seg = max(I_VE,I_seg);
    end
end

%% Remove segments with fewer than "thresh" voxels
function I_seg = clean_ves(I,thresh)

    CC = bwconncomp(I);
    I_seg = I;
    for uuu = 1:length(CC.PixelIdxList)
        if length(CC.PixelIdxList{uuu}) < thresh   % 30 for default
            I_seg(CC.PixelIdxList{uuu}) = 0;
        end
    end
end

%% cut away the side with only empty cell
function cellcube = cut_empty_cell(cellcube)
    cell_cube_notempty = cellfun(@isempty,cellcube)==0;
    idx = sum(cell_cube_notempty,[1,2])==0;
    cellcube(:,:,idx)=[];
    idx = sum(cell_cube_notempty,[2,3])==0;
    cellcube(idx,:,:)=[];
    idx = sum(cell_cube_notempty,[1,3])==0;
    cellcube(:,idx,:)=[];
end

%% Create matrix of logical True
function matrix_true = true_like(matrix)
    matrix_true = true(size(matrix));
end

%% Save MRI output
function save_mri_s(I, name, res, datatype)
    disp(' - making hdr...');
    % Make Nifti and header
    colres = res(2); 
    rowres = res(1); 
    sliceres = res(3); 
    % mri.vol = I;
    mri.volres = [res(1) res(2) res(3)];
    mri.xsize = rowres;
    mri.ysize = colres;
    mri.zsize = sliceres;
    a = diag([-colres rowres sliceres 1]);
    mri.vox2ras0 = a;
    mri.volsize = size(I);
    mri.vol = I;
    % mri.vol = flip(mri.vol,1);
    MRIwrite(mri,name,datatype);
    disp(' - done - ');
end

%% Calcualte dimensions for partitioning volume
function [xarray, yarray, zarray] = calc_part_dim(vol)
% Calculate partition dimensions
% Equally divide each dimension of the volume.
% TODO: ensure the partitions will fit into available memory.
%
% INPUTS:
%		vol (matrix): entire mri volume
% OUTPUTS:
%		pdims (matrix):  three arrays

%%% Find memory available for jobs
% Check available memory
m = memory;
% Convert bytes to gigabytes
m_gb = (m.MemAvailableAllArrays) / 1e9;

%%% Calculate partitions
% Find x,y,z dimensions of entire volume
[x,y,z] = size(vol);
% Calculate parition size
xarray = dimpart(x,0);
yarray = dimpart(y,0);
zarray = dimpart(z,1);

    %% Find the partition size of each dimension
    function [parray] = dimpart(dim, zbool)
    % Calculate the partition size for the respective dimension.
    %
    % INPUTS:
    %       dim (double): size of dimension
    %       zbool (double): boolean for z dimension
    % OUTPUTS
    %       parray (array): array containing an array for each dimension
    
    % If zbool is true, then dividend equals 2. This is because the z
    % dimension is always smaller than x,y, so it can safely be separated
    % into just two paritions.
    if zbool
        div = 2;
    else
        div = 3;
    end
    % Divide dimension by 3 to find partition dimension
    pdim = floor(dim ./ div);
    % Find remainder
    prem = rem(dim, div);
    % Create array for respective dimension of partition
    if zbool
        parray = [pdim, (pdim+prem)];
    else        
        parray = [pdim, pdim, (pdim+prem)];
    end
    
    end

end















