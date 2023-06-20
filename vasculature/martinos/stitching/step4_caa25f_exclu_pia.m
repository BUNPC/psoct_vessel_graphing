%% 

% Variables to pass in as input
%	

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

%% Create WM mask
nii = MRIread(mus_AIP_path);
figure;histogram(nii.vol(nii.vol~=0))
thres = multithresh(nii.vol,2);
mask = nii.vol > thres(2);
mask_filt = bwmorph3(mask,'majority');
mask = mask_filt;

%%% Select first two and last two slices
% TODO: programatically select first two slices & last two slices
mask(:,:,1) = imfill(mask(:,:,1),'holes');
mask(:,:,2) = imfill(mask(:,:,2),'holes');
mask(:,:,55) = imfill(mask(:,:,55),'holes');
mask(:,:,56) = imfill(mask(:,:,56),'holes');
% floating slices - set mask equivalent to nearest non-floating slice
mask(:,:,43) = imfill(mask(:,:,43),'holes');
mask(:,:,44) = imfill(mask(:,:,43),'holes'); % floating slice
mask(:,:,45) = imfill(mask(:,:,45),'holes'); % floating slice
mask(:,:,46) = imfill(mask(:,:,45),'holes');
% Completely fill hols
mask = imdilate(mask,strel('sphere',5));
mask = bwmorph3(mask,'fill');
mask = imfill(mask,'holes');
mask = imerode(mask,strel('sphere',5));
view3D(mask*1)
% Find largest closed object in mask -> convert to binary
cc = bwconncomp(mask);
vxsize=cellfun(@length,cc.PixelIdxList,'UniformOutput',true);
[~,idx] = max(vxsize);
mask_tmp = mask*0;
mask_tmp(cc.PixelIdxList{idx}) = 1;
view3D(mask_tmp)
mask = mask_tmp;

%% Create pia mask (Stacked_AIP_pia.nii.gz) (AIP conv threshold)
% Load AIP
nii_test = MRIread(AIP_path);
%%% Edge detection via convolution with 60x60 kernel.
% Method: - 2D convolution to each slice. 
%	      - tmp2 =  original slice divided by convolved slice.
tmp2 = conv3D(nii_test.vol,60);
% Account for floating slices
tmp2(:,:,44) = tmp2(:,:,43);
tmp2(:,:,45) = tmp2(:,:,46);
view3D(tmp2)

%%% take minimum intensity projection (minIP) of average intensity projection (AIP)
% ds is the number of slices to use for minIP
ds = 5;
minIP = do_Z_projection(ds,tmp2,'min');
view3D(minIP)
% tmp is the thresholded minIP
tmp = [minIP<0.96]*1;

%%% Smooth, threshold, close minIP. set WM = false.
tmp = imfilter(tmp,fspecial('gaussian',100,10));
mask_pia = tmp>0.2;
mask_pia = imclose(mask_pia,strel('sphere',5)); 
% set WH equal to false
mask_pia(mask) = false;

%%% Account for downsampling
mask_pia(:,:,end-1:end) = [];                                   	% remove last two slice
mask_pia = cat(3, mask_pia(:,:,1), mask_pia(:,:,1), mask_pia);     	% add two slice in front
view3D(mask_pia*1);
% TODO: pass in resolution
save_mri(mask_pia,'Stacked_AIP_pia.nii.gz',[0.01 0.01 0.1],'uchar')

%% Manually compare 
% Stack 1: convolved / normalized stack
% Stack 2: minimumIP of convolved/normalized output
% Ensure boundaries are enclosed around pia

%% [Freeview]  Manually refine the pia and tissue edge mask
% Note: this step may not be necessary. May be able to just clean agarose mask
%  input file = 'Stacked_AIP_pia.nii.gz'
%  output file = 'Stacked_AIP_pia.cleaned.nii.gz'

%% create agarose mask ( MIP -> agar -> nontissue )
MIP_path = ['/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/StackNii/'...
    'Stacked_MIP.nii'];
nii_test = MRIread(MIP_path);
view3D(nii_test.vol)
figure;
histogram(nii_test.vol(nii_test.vol~=0));
view3D(nii_test.vol>70);
mask_agar = imfill(nii_test.vol>70,'holes');

% Exclude floating slices
mask_agar(:,:,44) = mask_agar(:,:,43);
mask_agar(:,:,45) = mask_agar(:,:,46);
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
    'Stacked_AIP_pia.nii.gz']; nii_test = MRIread(pia_path);
mask_nontissue = [nii_test.vol | mask_agar];
% view3D(mask_nontissue)
save_mri(mask_nontissue,...
    '/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/StackNii/Stacked_nontissue.nii.gz',...
    [0.01 0.01 0.1],'uchar')

%% Create visual aid: apply tissue mask to Frangi filter output
% TODO: need method to determine which slices to exclude.
% Currently, Dylan asks Hui for guidance.
% Need exclusion criteria.
nii = MRIread('mus_mean_10um-iso.slice40px.I_seg.nii');
% Exclude floating slices (adjacent slices that received artifact from Gaussian smoothing)
nii.vol(:,:,431:450) = false;
% Exclude last two slices
nii.vol(:,:,561:end) = false;
% Maximum intensity projection (10 slices to match non-tissue mask)
I = do_downsample(nii.vol,[1,1,10],'max');
% TODO: pass in variable resolution.
save_mri(I,...
    '/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/StackNii/Stacked_seg_MIP.nii.gz',...
    [0.01 0.01 0.1],'uchar')
	
%% [Freeview]  Manually refine the nontissue mask (final revision)
% load "Stacked_nontissue.nii.gz" (this will be edited)
% overlay Stacked_seg_MIP.nii.gz and Stacked_AIP.nii
% output = Stacked_nontissue.cleaned.nii.gz
% This step requires manual intervention.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%% The following code is run after the manual cleaning
% This should be separate function.

% load cleaned nontissue mask
nii = MRIread('/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/StackNii/Stacked_nontissue_05302023.mgz');
mask_nontissue = nii.vol;
mask_nontissue = mask_nontissue ==1;

% Apply 10iso mask (pia + agarose) (repmat from [10/10/100] to 10iso)
nii_test = MRIread('/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/mus_vessel/mus_mean_10um-iso.slice40px.I_seg.nii');
I_seg = nii_test.vol ==1; 
mask_nontissue_rep = reshape(repmat(mask_nontissue,[1,10,1]),size(I_seg));
I_seg(mask_nontissue_rep)=false;

% Remove segmentation of slices w/ artifacts (+ adjacent slices)
I_seg(:,:,428:453)=false;
I_seg(:,:,541:end)=false;

% Delete temporary variable
clear mask_nontissue_rep
nii_test.vol = I_seg;
MRIwrite(nii_test,['/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/StackNii/' ... 
    'mus_mean_10um-iso.slice40px.I_seg.excl_pia.nii'],'uchar');
% Delete temporary variable
clear I_seg

% Create z projection for QC
% imwrite(max(nii_test.vol(:,:,1:400),[],3),'ves_MIP.before.z400.tiff')
% imwrite(max(I_seg(:,:,1:400),[],3),'ves_MIP.after.z400.tiff')

%% Create graph.mat 
% TODO: change this to run function.
%  from 'mus_mean_10um-iso.slice40px.I_seg.excl_pia.clean.nii'
run('/autofs/cluster/octdata2/users/Hui/tools/dg_utils/ves2graph/SegtoGraph.m');
save(['/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/mus_vessel/' ...
    'mus_mean_10um-iso.slice40px.I_seg.excl_pia.cleaned.graph.mat'],'Graph');

%% save skeleton as NIFTI
load('/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/mus_vessel/mus_mean_10um-iso.slice40px.I_seg.excl_pia.cleaned.graph.mat')
name = '/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/mus_vessel/mus_mean_10um-iso.slice40px.I_seg.excl_pia.cleaned.graph.sk.nii';
sk3D([2879, 2613, 690],Graph,name,[0.01,0.01,0.01],0,1);

%% Move these functions to separate function files
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