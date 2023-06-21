
%% load
mri = MRIread(['/autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/process_run1/mus_vessel' ...
               '/mus_mean_10um-iso.slice40px.nii']);
fprintf('mri loaded.\n')
mri.vol = single(mri.vol);

%% divide volume
cell_I = mat2cell(mri.vol,[1000,1000,1099],[1000,1000,1410],[600,600]);


%%
% padded cell variable
cell_I_padded = cell(size(cell_I)+2);
cell_I_padded = cellfun(@single,cell_I_padded,'UniformOutput',false);
cell_I_padded2 = cell_I_padded;
cell_I_m_padded = cell_I_padded;
cell_I_padded(2:end-1,2:end-1,2:end-1) = cell_I; % save partitioned images in the center


% set skirt range str [20] 
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



%% save partition (maybe remove?)
for s = 1:length(cell_I2(:))
    var_name = sprintf('I%i',s);
    eval([var_name ' = cell_I2{s};'])
    save('Ma_partition_cell_padded.mat',var_name,'-append')
    eval(['clear ' var_name])
end
reshape_sz = size(cell_I2);
save('Ma_partition_cell_padded.mat','reshape_sz','-append')
save('Ma_partition_cell_padded.mat','cell_I2_skirt_lnridx','-append')
save('Ma_partition_cell_padded.mat','cell_I2_main','-append')
save('Ma_partition_cell_padded.mat','cell_I2_main_skirt','-append')


%% Image Processing and Frangi Filter
% TODO: call my Frangi filter function
pause(0.1)

%% I_seg reconstruct
load('Ma_partition_cell_padded.mat','reshape_sz','cell_I2_main',...
    'cell_I2_main_skirt','cell_I2_skirt_lnridx');
cell_seg = cell(reshape_sz);

for s = 1:length(cell_I2_main(:))
    var_in_name = sprintf('I_seg%i',s);
    load('ves4_padded.mat',var_in_name)
    eval(['I_seg = ' var_in_name ';'])
    eval(['clear ' var_in_name])

    I_seg(cell_I2_skirt_lnridx{s})=[];
    I_seg = reshape(I_seg, cell_I2_main{s});
    cell_seg{s} = I_seg;
end


I_seg_all = cell2mat(cell_seg);
save_mri(I_seg_all, 'mus_mean_10um-iso.slice40px.I_seg.nii', [0.01,0.01,0.01],'uchar')

%% create_mri_mosaic.py commandline 
%  create_mri_mosaic.py is Jackson's script to divide large nii image into 
%  small blocks using freesurfer functions. The small blocks retained its
%  *original positions* when loaded in freeview. 
%% create_mri_mosaic.py -> I_seg
% /autofs/space/tiamat_001/users/jackson/projects/scripts/Python/create_mri_mosaic.py \ 
% /autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/process_run1/mus_vessel/mus_mean_10um-iso.slice40px.I_seg.nii 1100 1000 1200 --m_dims 3 3 1   -s --outdir /autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/process_run1/mus_vessel/nii_split --tile_name seg
%% create_mri_mosaic.py -> I
% /autofs/space/tiamat_001/users/jackson/projects/scripts/Python/create_mri_mosaic.py \
% /autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/process_run1/mus_vessel/mus_mean_10um-iso.slice40px.t3_9.sigma1.nii 1100 1000 1200 --m_dims 3 3 1   -s --outdir /autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/process_run1/mus_vessel/nii_split --tile_name I

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

function matrix_true = true_like(matrix)
    matrix_true = true(size(matrix));
end

