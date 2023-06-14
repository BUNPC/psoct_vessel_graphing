% mus_minIP5_path = ['/autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/process_run1/StackNii/' ...
%     'Stacked_mus_minIP5.nii'];
mus_AIP_path = ['/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/StackNii/' ...
    'Stacked_mus_AIP.nii'];
% mus_minIP5_conv3_path = ['/autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/process_run1/StackNii/'...
%     'Stacked_minIP1.conv3.minIP5.nii']; nii_test = MRIread(mus_minIP5_conv3_path);
% mus_minIP_conv3_path = ['/autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/process_run1/StackNii/'...
%     'Stacked_minIP1.conv3.nii']; nii_test = MRIread(mus_minIP_conv3_path);
AIP_path = ['/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/StackNii/'...
    'Stacked_AIP.nii']; 
%%
%% Create WM mask
nii = MRIread(mus_AIP_path);
figure;histogram(nii.vol(nii.vol~=0))
thres = multithresh(nii.vol,2);
mask = nii.vol > thres(2);

mask_filt = bwmorph3(mask,'majority');
mask = mask_filt;

mask(:,:,1) = imfill(mask(:,:,1),'holes');
mask(:,:,2) = imfill(mask(:,:,2),'holes');

mask(:,:,43) = imfill(mask(:,:,end),'holes');
mask(:,:,44) = imfill(mask(:,:,end),'holes'); % floating slice
mask(:,:,45) = imfill(mask(:,:,end),'holes'); % floating slice
mask(:,:,46) = imfill(mask(:,:,end),'holes');
mask(:,:,55) = imfill(mask(:,:,end),'holes');
mask(:,:,56) = imfill(mask(:,:,end),'holes');

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
tmp2(:,:,44) = tmp2(:,:,43); % floating slice
tmp2(:,:,45) = tmp2(:,:,46); % floating slice
view3D(tmp2)
ds = 5;
minIP2 = do_Z_projection(ds,tmp2,'min');
view3D(minIP2)

tmp = [minIP2<0.96]*1;
% for repeat = 1:4;tmp = imfilter(tmp,fspecial('gaussian',30,10));end
tmp = imfilter(tmp,fspecial('gaussian',100,10));
mask_pia = tmp>0.2;
mask_pia = imclose(mask_pia,strel('sphere',5)); 
mask_pia(mask) = false; % mask out anything in the WM
mask_pia(:,:,end-1:end) = [];                                   % remove last two slice
mask_pia = cat(3,mask_pia(:,:,1),mask_pia(:,:,1),mask_pia);     % add two slice in front
view3D(mask_pia*1);%figure; imshow(tmp);
save_mri(mask_pia,'Stacked_AIP_pia.nii.gz',[0.01 0.01 0.1],'uchar')
%% [Freeview]  Manually refine the pia and tissue edge mask
%  output is named 'Stacked_AIP_pia.cleaned.nii.gz'
%% create agarose mask ( MIP -> agar -> nontissue )
MIP_path = ['/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/StackNii/'...
    'Stacked_MIP.nii']; nii_test = MRIread(MIP_path);
view3D(nii_test.vol)
figure; histogram(nii_test.vol(nii_test.vol~=0))
view3D(nii_test.vol>70)
mask_agar = imfill(nii_test.vol>70,'holes');

mask_agar(:,:,44) = mask_agar(:,:,43); % floating slice
mask_agar(:,:,45) = mask_agar(:,:,46); % floating slice
mask_agar = get_largest(mask_agar);
% mask_agar = bwmorph3(mask_agar,'majority');
mask_agar = imdilate(mask_agar,strel('sphere',5));
mask_agar = bwmorph3(mask_agar,'fill');
mask_agar = imerode(mask_agar,strel('sphere',5));
mask_agar = ~mask_agar;
mask_agar = get_largest(mask_agar);
view3D(mask_agar)
mask_agar(:,:,44:45) = true;    % consider slice 44 and 45 as nontissue
mask_agar(:,:,57:end) = true;   % consider slice 57 and after as nontissue
%% Create Stacked_nontissue.nii.gz
pia_path = ['/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/StackNii/'...
    'Stacked_AIP_pia.cleaned.nii.gz']; nii_test = MRIread(pia_path);
mask_nontissue = [nii_test.vol | mask_agar];
% view3D(mask_nontissue)
save_mri(mask_nontissue,...
    '/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/StackNii/Stacked_nontissue.nii.gz',...
    [0.01 0.01 0.1],'uchar')

%%
nii = MRIread('mus_mean_10um-iso.slice40px.I_seg.nii');
nii.vol(:,:,431:450) = false;
nii.vol(:,:,561:end) = false;
I = do_downsample(nii.vol,[1,1,10],'max');
save_mri(I,...
    '/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/StackNii/Stacked_seg_MIP.nii.gz',...
    [0.01 0.01 0.1],'uchar')
%% load cleaned nontissue mask
nii = MRIread('/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/StackNii/Stacked_nontissue_05302023.mgz');
mask_nontissue = nii.vol;
mask_nontissue = mask_nontissue ==1;

%% Apply 10iso mask (pia + agarose) (repmat from [10/10/100] to 10iso)
nii_test = MRIread('/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/mus_vessel/mus_mean_10um-iso.slice40px.I_seg.nii');
I_seg = nii_test.vol ==1; 
mask_nontissue_rep = reshape(repmat(mask_nontissue,[1,10,1]),size(I_seg));
I_seg(mask_nontissue_rep)=false;
I_seg(:,:,428:453)=false; % slice 44 and 45 are zeroed, adjacent planes are zeroed as well
I_seg(:,:,541:end)=false; % slice 55 and onward are zeroed
clear mask_nontissue_rep
nii_test.vol = I_seg;
MRIwrite(nii_test,['/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/StackNii/' ... 
    'mus_mean_10um-iso.slice40px.I_seg.excl_pia.nii'],'uchar');
clear I_seg

% imwrite(max(nii_test.vol(:,:,1:400),[],3),'ves_MIP.before.z400.tiff')
% imwrite(max(I_seg(:,:,1:400),[],3),'ves_MIP.after.z400.tiff')
%% [Freeview]  Manually clean-up floating slice artifact and tissue edge 
%  output is named 'mus_mean_10um-iso.slice40px.I_seg.excl_pia.clean.nii'
%% Create graph.mat 
%  from 'mus_mean_10um-iso.slice40px.I_seg.excl_pia.clean.nii'
run('/autofs/cluster/octdata2/users/Hui/tools/dg_utils/ves2graph/SegtoGraph.m');
save(['/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/mus_vessel/' ...
    'mus_mean_10um-iso.slice40px.I_seg.excl_pia.cleaned.graph.mat'],'Graph');


%% save sk for whole tissue no ds
load('/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/mus_vessel/mus_mean_10um-iso.slice40px.I_seg.excl_pia.cleaned.graph.mat')
name = '/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/mus_vessel/mus_mean_10um-iso.slice40px.I_seg.excl_pia.cleaned.graph.sk.nii';
sk3D([2879, 2613, 690],Graph,name,[0.01,0.01,0.01],0,1);


%%
% nii_test = MRIread('/autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/process_run1/mus_vessel/mus_mean_10um-iso.slice40px.t3_9.sigma1.nii');
% MAT2TIFFuint8(1-nii_test.vol, '/autofs/cluster/octdata2/users/Chao/caa/caa_17/occipital/process_run1/mus_vessel/mus_mean_10um-iso.slice40px.t3_9.sigma1.tiff')
% tiff format has size limit and our tissue has exceeded that. 
vesGraphValidate
%%
function I_seg = get_largest(I_seg)
    CC = bwconncomp(I_seg,18); I_seg(:) = false; 
    [~,idx] = max(cellfun(@length,CC.PixelIdxList));
    I_seg(CC.PixelIdxList{idx}) = true;
end
function I_seg = clean_ves(I,thresh)

    CC = bwconncomp(I);
    I_seg = I;
    for uuu = 1:length(CC.PixelIdxList)
        if length(CC.PixelIdxList{uuu}) < thresh   % 30 for default
            I_seg(CC.PixelIdxList{uuu}) = 0;
        end
    end
end