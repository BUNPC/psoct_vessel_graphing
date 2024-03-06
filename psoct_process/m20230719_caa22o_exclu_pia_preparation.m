addpath('/autofs/cluster/octdata2/users/Hui/tools/dg_utils/vol_recon_beta');
% mus_minIP5_path = ['/autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/process_run1/StackNii/' ...
%     'Stacked_mus_minIP5.nii'];
mus_AIP_path = ['/autofs/cluster/octdata2/users/Chao/caa/caa_22/processed/20211007/StackNii/' ...
    'Stacked_mus_AIP.nii'];
% mus_minIP5_conv3_path = ['/autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/process_run1/StackNii/'...
%     'Stacked_minIP1.conv3.minIP5.nii']; nii_test = MRIread(mus_minIP5_conv3_path);
% mus_minIP_conv3_path = ['/autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/process_run1/StackNii/'...
%     'Stacked_minIP1.conv3.nii']; nii_test = MRIread(mus_minIP_conv3_path);
AIP_path = ['/autofs/cluster/octdata2/users/Chao/caa/caa_22/processed/20211007/StackNii/'...
    'Stacked_AIP.nii']; 

MIP_path = ['/autofs/cluster/octdata2/users/Chao/caa/caa_22/processed/20211007/StackNii/'...
    'Stacked_MIP.nii']; 
%% generate Stacked_seg_MIP.nii.gz
nii = MRIread('/autofs/cluster/octdata2/users/Chao/caa/caa_22/processed/20211007/mus_vessel/mus_mean_10um-iso.slice40px.I_seg.nii.gz');
% nii.vol(:,:,431:450) = false;
% nii.vol(:,:,561:end) = false;
I = do_downsample(nii.vol,[1,1,10],'max');
%%
save_mri(I,...
    '/autofs/cluster/octdata2/users/Chao/caa/caa_22/processed/20211007/StackNii/Stacked_seg_MIP.DG2.nii.gz',...
    [0.01 0.01 0.1],'uchar')


%% Create WM mask
nii = MRIread(mus_AIP_path);
figure;histogram(nii.vol(nii.vol~=0))
thres = multithresh(nii.vol,2);
mask = nii.vol > thres(2);

mask_filt = bwmorph3(mask,'majority');
mask = mask_filt;

% mask(:,:,end) = imfill(mask(:,:,end),'holes');
mask = imdilate(mask,strel('sphere',5));
mask = bwmorph3(mask,'fill');
mask = imfill(mask,'holes');
mask = imerode(mask,strel('sphere',5));
view3D(mask*1)

cc = bwconncomp(mask);
vxsize=cellfun(@length,cc.PixelIdxList,'UniformOutput',true);
[~,idx] = max(vxsize);
mask_tmp = mask*0;
mask_tmp(cc.PixelIdxList{idx}) = 1;
view3D(mask_tmp)
mask = mask_tmp;

%% Create Stacked_AIP_pia.nii.gz (pia mask)( AIP conv threshold )
nii_test = MRIread(AIP_path);
tmp2 = conv3D(nii_test.vol,60);
view3D(tmp2)
ds = 5;
minIP2 = do_Z_projection(ds,tmp2,'min');
view3D(minIP2)

tmp = [minIP2<0.96]*1;
% for repeat = 1:4;tmp = imfilter(tmp,fspecial('gaussian',30,10));end
tmp = imfilter(tmp,fspecial('gaussian',100,10));
mask_pia = tmp>0.2;
mask_pia = imclose(mask_pia,strel('sphere',5)); 
% mask_pia(mask) = false; % mask out anything in the WM
mask_pia(:,:,end-1:end) = [];                                   % remove last two slice
mask_pia = cat(3,mask_pia(:,:,1),mask_pia(:,:,1),mask_pia);     % add two slice in front
view3D(mask_pia*1);%figure; imshow(tmp);
save_mri(mask_pia,'Stacked_AIP_pia.nii.gz',[0.01 0.01 0.1],'uchar')

%% create agarose mask ( MIP -> agar -> nontissue )
nii_test = MRIread(MIP_path);
view3D(nii_test.vol)
figure; histogram(nii_test.vol(nii_test.vol~=0))
view3D((nii_test.vol>61)*1)
mask_agar = imfill(nii_test.vol>61,'holes');
%%
mask_agar = get_largest(mask_agar);
% mask_agar = bwmorph3(mask_agar,'majority');
mask_agar = imdilate(mask_agar,strel('sphere',5));
mask_agar = bwmorph3(mask_agar,'fill');
mask_agar = imerode(mask_agar,strel('sphere',5));
mask_agar = ~mask_agar;
mask_agar = get_largest(mask_agar);
view3D(mask_agar)

%% Create Stacked_nontissue.nii.gz
pia_path = ['/autofs/cluster/octdata2/users/Chao/caa/caa_22/processed/20211007/StackNii/'...
    'Stacked_AIP_pia.nii.gz']; nii_test = MRIread(pia_path);
mask_nontissue = [nii_test.vol==1 | mask_agar];
% view3D(mask_nontissue)
save_mri(mask_nontissue,...
    '/autofs/cluster/octdata2/users/Chao/caa/caa_22/processed/20211007/StackNii/Stacked_nontissue.nii.gz',...
    [0.01 0.01 0.1],'uchar')

%% load cleaned nontissue mask
nii = MRIread('/autofs/cluster/octdata2/users/Chao/caa/caa_22/processed/20211007/StackNii/Stacked_nontissue.DJE.mgz');
mask_nontissue = nii.vol;
mask_nontissue = mask_nontissue ==1;

%% Apply 10iso mask (pia + agarose) (repmat from [10/10/100] to 10iso)
nii_test = MRIread('/autofs/cluster/octdata2/users/Chao/caa/caa_22/processed/20211007/mus_vessel/mus_mean_10um-iso.slice40px.I_seg.nii.gz');
I_seg = nii_test.vol ==1; 
mask_nontissue_rep = reshape(repmat(mask_nontissue,[1,10,1]),size(I_seg,1),size(I_seg,2),[]);

I_seg(mask_nontissue_rep)=false;
clear mask_nontissue_rep
nii_test.vol = I_seg;
MRIwrite(nii_test,['/autofs/cluster/octdata2/users/Chao/caa/caa_22/processed/20211007/mus_vessel/' ... 
    'mus_mean_10um-iso.slice40px.I_seg.excl_pia.mgz'],'uchar');
clear I_seg

imwrite(max(nii_test.vol(:,:,1:400),[],3),'/autofs/cluster/octdata2/users/Chao/caa/caa_22/processed/20211007/mus_vessel/ves_MIP.before.z400.tiff')
% imwrite(max(I_seg(:,:,1:400),[],3),'ves_MIP.after.z400.tiff')
%% Create graph.mat 
addpath(genpath('/autofs/cluster/octdata2/users/Hui/tools/dg_utils/ves2graph'))
%  from 'mus_mean_10um-iso.slice40px.I_seg.excl_pia.clean.nii'
run('/autofs/cluster/octdata2/users/Hui/tools/dg_utils/ves2graph/SegtoGraph.m');
save(['/autofs/cluster/octdata2/users/Chao/caa/caa_22/processed/20211007/mus_vessel/' ...
    'mus_mean_10um-iso.slice40px.I_seg.excl_pia.cleaned.graph.mat'],'Graph');


%% save sk for whole tissue no ds
load(['/autofs/cluster/octdata2/users/Chao/caa/caa_22/processed/20211007/mus_vessel/' ...
    'mus_mean_10um-iso.slice40px.I_seg.excl_pia.cleaned.graph.mat'],'Graph');
name = '/autofs/cluster/octdata2/users/Chao/caa/caa_22/processed/20211007/mus_vessel/mus_mean_10um-iso.slice40px.I_seg.excl_pia.cleaned.graph.sk.mgz';
sk3D([2561,3926,700],Graph,name,[0.01,0.01,0.01],0,1);

%%
function I_seg = get_largest(I_seg)
    CC = bwconncomp(I_seg,18); I_seg(:) = false; 
    [~,idx] = max(cellfun(@length,CC.PixelIdxList));
    I_seg(CC.PixelIdxList{idx}) = true;
end