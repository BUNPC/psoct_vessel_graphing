%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stitching Z
load(ParameterFile)

sliceidx = Mosaic3D.sliceidx;
modality = 'dBI';
zoff        = Mosaic3D.z_parameters(1);
z_sm        = Mosaic3D.z_parameters(2);
z_s         = Mosaic3D.z_parameters(3);
z_sms       = z_sm+z_s; 
z_m         = z_sm-z_s;

fprintf('zoff    skirt      mainbody    skirt\n')
fprintf('%0.2i      %0.2i         %0.2i          %0.2i\n',zoff, z_s,z_m,z_s)
fprintf('%0.2i + %0.2i = %0.2ipx (%ium) -> slice thickness\n',z_s,z_m,z_sm,z_sm*2.5)

%zoff    skirt      mainbody    skirt
%5       15         25          15

%skirt + mainbody = slice thickness
%15 + 25 = 40px (100um)

load([Scan.ProcessDir filesep modality filesep modality '_slice003_mean-xy_20um-iso.mat'],'Id');%c

if strcmpi(modality,'mus')
    mus_aip = mean(smooth3(Id),3);
    thresh = multithresh(mus_aip,2); % output will be [agar-gm gm-wm]
    figure;subplot(1,2,1);imshow(mus_aip,[thresh(2),prctile(mus_aip(:),98)]); title('wm');
    subplot(1,2,2);imshow(mus_aip,[thresh(1),thresh(2)]); title('gm');


    wm_mask = mus_aip>thresh(2);    mask3d = logical(Id(1,1,:)*0 + wm_mask);
    % gm_mask = mus_aip<thresh(2) & mus_aip>thresh(1);    mask3d = logical(Id(1,1,:)*0 + gm_mask);

    tissue=mean(reshape(Id(mask3d),[],size(Id,3)), [1]);clear mask3d
    tissue=tissue(:).';
    figure; plot(tissue)
elseif strcmpi(modality,'dBI')
    % choose a tissue-only block to generate 3 weights
    tissue=mean(Id(850:850+190,1370:1370+400,1+zoff:end), [1 2]);
    tissue=tissue(:).';
end
% choose a tissue-only block to generate 3 weights
%tissue=squeeze(mean(Id(100:100+300,600:600+300,1+zoff:end), [1 2])).'; %caa6 occi
% tissue=squeeze(mean(Id(850:850+190,1370:1370+400,1+zoff:end), [1 2])).'; %view3D(Id(250:250+300,870:600+300,1+zoff:end))

% generate weights that uniform intensity level across depth
tissue = tissue([1:z_sms]+zoff);
s  = tissue(z_s+1)./tissue(1:z_s);        % Top Overlapping Skirt
ms = tissue(z_s+1)./tissue(z_s+1:z_sms);  % Non-overlapping Body  &  Bottom Overlapping Skirt

% blending weight at the overlap
switch modality
    case 'dBI'; degree = 1;   
    case 'mus'; degree = 1; % 1  
end
w1 = s             .*linspace(0,1,z_s).^degree; w1 = w1(:).';    % Top Overlapping Skirt
w2 = ms(z_m+1:end) .*linspace(1,0,z_s).^degree; w2 = w2(:).';    % Bottom Overlapping Skirt
w3 = ms(1:z_m);                                 w3 = w3(:).';


nb_slices   = size(sliceidx,2);
[row_pxblock, col_pxblock] = size(Id,[1 2]);
tot_z_pix   = z_sm*nb_slices;

    %%% template for stitched block
Ma = zeros(row_pxblock, col_pxblock, tot_z_pix,'single');
%     Ma = zeros(w, l, onl*(nb_slices));

clear Id s ms
%% Begin stitching  
tic

indir = [Scan.ProcessDir filesep modality];
outdir = indir;
fprintf(' - Input directory = \n%s\n',indir);
fprintf(' - Output directory = \n%s\n',outdir);

% w1=
% w3=
% w2=

 
    
%%% loop slices
for s = 1:size(sliceidx,2)

    
    %kk_real = kk-first_slice+1; %kk+1;
    fprintf('\tslice %d\n',s);

    %%% load slice
    fdata = [indir filesep modality '_slice' sprintf('%03i',s) '_mean-xy_20um-iso.mat'];
    load(fdata,'Id');

    %eval(['data' num2str(kk) '=data;']);
    ddims = size(Id);

    %%% stitch (cases for first, last, and other slice #s)
    if s==1 % first_slice
        W = reshape(repmat([w1*0+1 w3 w2],          [ddims(1)*ddims(2),1]),ddims(1),ddims(2),[]);
        zr1 = 1:z_sms;
        zr2 = zr1+zoff; %zr2 = zr1; %
        Ma(:,:,zr1) = Id(:,:,zr2).*W;

    elseif s==nb_slices % last_slice
        % Add to the `t` bottom-most pixels in Ma
        W = reshape(repmat([w1 w3 ],          [ddims(1)*ddims(2),1]),ddims(1),ddims(2),[]);            
        zr1 = [1:z_sm] + tot_z_pix-z_sm;
        zr2 = zoff+[1:z_sm]; %zr2 = [1:z_sms]; %
        Ma(:,:,zr1) = Ma(:,:,zr1) + Id(:,:,zr2).*W;
    else
        % For slices (first_slice+1) to (last_slice-1)
        W = reshape(repmat([w1 w3 w2],              [ddims(1)*ddims(2),1]),ddims(1),ddims(2),[]);
        zr1 = [1:z_sms] + (s-1)*z_sm;
        zr2 = zoff+[1:z_sms]; %zr2 = [1:z_sms]; %
        Ma(:,:,zr1) = Ma(:,:,zr1) + Id(:,:,zr2).*W;
    end
%         disp('Finished slice');

%         subplot(10,1,kk);
%         imshow(squeeze(Ma(:,323,:)),[0 0.25]); colormap jet; drawnow;

end
%     fprintf('Saving block %d stitching...\n',currblock);
%     
%     
%     name = [outdir filesep modality '_stitch_block' sprintf('%03i',currblock) '.nii'];
%     resolution = [in_planeRes,in_planeRes,thru_planeRes];
% %     MakeNiiLAS(name,Ma,resolution);
%     
%     disp(' - making hdr...');
%     % Make Nifti and header
%     colres = resolution(2); 
%     rowres = resolution(1); 
%     sliceres = resolution(3); 
%     % mri.vol = I;
%     mri.volres = [resolution(1) resolution(2) resolution(3)];
%     mri.xsize = rowres;
%     mri.ysize = colres;
%     mri.zsize = sliceres;
%     a = diag([-colres rowres sliceres 1]);
%     mri.vox2ras0 = a;
%     mri.volsize = size(Ma);
%     mri.vol = Ma;
%     % mri.vol = flip(mri.vol,1);
%     MRIwrite(mri,name,'float');


toc

    %         view3D(Id);

%%
    
    disp(' - Downsampling to 20 micron iso (z)...');
    %m=1:sxy:size(Ma,1)-sxy;
    %n=1:sxy:size(Ma,2)-sxy;
    p=1:8:size(Ma,3)-8+1; %size(I,3);
    %Id = zeros(length(m),length(n),length(p));
    
%     sx = 1; %20um
%     sy = 1; %20um
%     sz = 8; %2.5um -> 20um
    ds_xyz = [1, 1, 8]; % sz are loaded from recon.mat

    [Id] = do_downsample(Ma, ds_xyz);
    
    resolution      = [0.01 0.01 0.0025];
    resolution_ds = resolution .* [2 2 8];
    
    name = [outdir filesep modality '_mean_20um-iso.nii'];

    disp(' - making hdr...');
    % Make Nifti and header
    colres = resolution_ds(2); 
    rowres = resolution_ds(1); 
    sliceres = resolution_ds(3); 
    % mri.vol = I;
    mri.volres = [resolution_ds(1) resolution_ds(2) resolution_ds(3)];
    mri.xsize = rowres;
    mri.ysize = colres;
    mri.zsize = sliceres;
    a = diag([-colres rowres sliceres 1]);
    mri.vox2ras0 = a;
    mri.volsize = size(Id);
    mri.vol = Id;
    % mri.vol = flip(mri.vol,1);
    MRIwrite(mri,name,'float');
    
