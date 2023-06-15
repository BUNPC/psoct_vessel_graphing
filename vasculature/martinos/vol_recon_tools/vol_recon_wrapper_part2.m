clear 

cd /autofs/cluster/octdata2/users/Hui/tools/dg_utils/vol_recon_beta
%% Mosaic 3D Slices
modality = 'mus';                                                           % 'mus','dBI';

outdir = [para.ProcDir{1} filesep modality];
MZL = 30;
sxyz = [2, 2, 8]; 
tic
%poolobj = parpool(5);
for sliceid     = para.out %[19:21] %sliceindex(2,:)
    sliceinfo   = sliceidx(:,sliceid);
    rawdatadir  = para.rawdatadir{sliceinfo(3,1)};
    
    Mosaic_Telesto_3D_Fiji_script(rawdatadir, outdir, expFile, modality, sliceinfo, MZL, sxyz);
end
%delete(poolobj)
toc


%% Recon Parameter file
%z_stitching parameter
first_slice = 1;%slice_range{run_n}(1);
last_slice = 120;
nb_slices = last_slice;

z_sms = 25;         % 's' refers to the skirt of the slice head/tail. 'm' refers to main body of the slice
z_sm = 20;          % slice thickness(px)
z_s = z_sms - z_sm;           % z_sms - slice thickness(px) 
z_m = z_sms - 2*z_s;
tot_z_pix = z_sms + (nb_slices-1)*z_sm;  %%% total z pixels in stitched block

zoff = 5;    %pixel offset in z depth (surface crop)
nz = 25;
z_range = [zoff+1]:[zoff+nz];


% row (dim 1) & col (dim 2) (og: 3117 x 3410)  (ds: 1558 x 1705)
row_s           = 1;                        % dim=1 in matlab
row_e           = 1558;
row_range       = row_s:row_e;
row_pixels_tot  = length(row_range);
row_nblock      = 2; 
row_pxblock     = row_pixels_tot/row_nblock;   

col_s           = 1;                        % dim=2 in matlab  
col_e           = 1705;
col_range       = col_s:col_e;
col_pixels_tot  = length(col_range);
col_nblock      = 5;          
col_pxblock     = col_pixels_tot/col_nblock;      


% x_sz * y_sz * sm * nSlice * 4bit = rough file size (double)

% Resolution & Dimension

%%% FOR orignal
    thru_planeRes   = 0.0025; % it was 0.0025, but i think it should be 0.0029... % z resolution
    in_planeRes     = 0.01; % xy-resolution
    fprintf(' xy res=%g mm\n z  res=%g mm\n',in_planeRes,thru_planeRes);
    resolution      = [in_planeRes in_planeRes thru_planeRes];
% 
%     origdims        = [row_pixels_tot col_pixels_tot tot_z_pix];


%%% FOR downsampling (ds)
    sxy = 2;
    sz  = 8;
    resolution_ds = resolution .* [sxy sxy sz];

%     m=1:sxy:row_pixels_tot;
%     n=1:sxy:col_pixels_tot;
    p=1:sz:tot_z_pix-sz+1; %size(I,3);
    dsdims   = [length(row_range),length(col_range),length(p)]; %[length(m) length(n) length(p)];

    row_pxblock_ds = row_pxblock; %row_pxblock/sxy; %dim=1
    col_pxblock_ds = col_pxblock; %col_pxblock/sxy; %dim=2

save(reconFile, ...
    'first_slice', 'last_slice', 'nb_slices',...
    'z_sms', 'z_s', 'z_sm', 'z_m', 'tot_z_pix', 'zoff', 'z_range',... 
    'row_s', 'row_e', 'row_range', 'row_pixels_tot', 'row_nblock', 'row_pxblock',... %
    'col_s', 'col_e', 'col_range', 'col_pixels_tot', 'col_nblock', 'col_pxblock',... % 
    'thru_planeRes', 'in_planeRes', 'resolution',... 'origdims',
    'sxy', 'sz', 'resolution_ds', 'dsdims',  'row_pxblock_ds', 'col_pxblock_ds');... 'm', 'n','p', 
     %'ProcDir'

%% BlockDivision
ProcDir=para.ProcDir{1};
modality = 'dBI';           % 'mus','dBI';

moddir = [ProcDir filesep modality];
indir = moddir; %[moddir filesep '1_xy_Mosiac']
outdir = [moddir filesep 'XY_Division'];

tic
poolobj = parpool(5);
parfor sliceid = 1:120 
    sliceid

    do_xy_Division_ds( indir, outdir, reconFile, modality, sliceid);
    
end
toc
delete(poolobj)  

%% BlockStitching
modality = 'dBI';           % 'mus','dBI';
ProcDir=para.ProcDir{1};
% choose a tissue-only block to generate 3 weights
fdata = [ProcDir filesep modality filesep modality '_slice003_mean-xy_20um-iso.mat'];
load(fdata,'Id');%view3D(data);
tissue=squeeze(mean(Id(600:end-600,600:end-600,:), [1 2])).'; clear Id

% generate weights that uniform intensity level across depth
s  = tissue(z_s)./tissue(1:z_s);        % Top Overlapping Skirt
ms = tissue(z_s)./tissue(z_s+1:z_sms);  % Non-overlapping Body  &  Bottom Overlapping Skirt

% blending weight at the overlap
switch modality
    case 'dBI'; degree = 1;   
    case 'mus'; degree = 1.8;   
end
w1 = s             .*linspace(0,1,z_s).^degree;     % Top Overlapping Skirt
w2 = ms(z_m+1:end) .*linspace(1,0,z_s).^degree;     % Bottom Overlapping Skirt
w3 = ms(1:z_m);                                     % Non-overlapping Body 


if 1==0
load(fdata,'Id'); view3D(Id); agar=squeeze(mean(Id(950:1050,1550:1650,:), [1 2])).'; clear data
z_stitching_display({tissue}, {'white matter'}, reconFile,w1,w2,w3);
end 

% blockstitching
moddir = [ProcDir filesep modality];
indir = [moddir filesep 'XY_Division'];
outdir = [ProcDir filesep 'Z_Stitched'];

poolobj = parpool(5);
parfor blockid = 1:10 %row_nblock*col_nblock
    do_z_Stitching(ProcDir, indir, outdir, reconFile, modality, blockid, w1,w2,w3)
end
delete(poolobj)

%% Merge
modality = 'dBI';           % 'mus','dBI';

indir = [ProcDir filesep 'Z_Stitched']; %[moddir filesep '1_xy_Mosiac']
outdir = indir;

MergeStitchedBlocks( indir, outdir, reconFile, modality )
