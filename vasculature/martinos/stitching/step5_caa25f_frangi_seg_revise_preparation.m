addpath(genpath('/autofs/cluster/octdata2/users/Hui/tools/dg_utils/ves2graph'))
addpath('/autofs/cluster/octdata2/users/Hui/tools/dg_utils/vol_recon_beta');
addpath('/autofs/cluster/octdata2/users/Hui/tools/dg_utils/matlab_fun')
%%
%  non-downsampled Graph is outputed
% 
%
%%
nii = MRIread('/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/mus_vessel/mus_mean_10um-iso.slice40px.nii');
mus = nii.vol;
mus(:,:,431:450) = 10; % to counter artifact [floating slice]
mus(:,:,571:end) = 10; % to ignore all data in this range [end of tissue]
mus = smooth3(mus,'gaussian',5,1.5);
view3D(mus)

nii.vol = mus;
MRIwrite(nii,'/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/mus_vessel/mus_mean_10um-iso.slice40px.sigma1p5.nii','float')

%% select tile from mosaic
[y,x] = meshgrid(1:3,1:3);
z = ones(size(y));
coords = [y(:),x(:),z(:)] - 1;

tilenum = 8;
tilecoords = coords(tilenum,:);

nii_dir = '/autofs/space/omega_001/users/caa/CAA25_Frontal/Process_caa25_frontal/mus_vessel/nii_split';
pad_coord = sprintf('mosaic_%i_%i_%i',tilecoords(:));
nii_name1 = sprintf('%s/seg_%s.mgz',nii_dir,pad_coord);
nii_name2 = sprintf('%s/I_%s.mgz',nii_dir,pad_coord);
nii_name3 = sprintf('%s/sk_%s.mgz',nii_dir,pad_coord);
nii_name4 = sprintf('%s/VE_%s.mgz',nii_dir,pad_coord);
%% run seg2graph
run('/autofs/cluster/octdata2/users/Hui/tools/dg_utils/ves2graph/SegtoGraph.m');
save([nii_dir '/graph_' pad_coord '.mat'],'Graph');

%% sk.nii.gz
name_sk = sprintf('%s/seg_%s.sk.nii.gz',nii_dir,pad_coord);
res = [0.01,0.01,0.01];
flag_downsampled = 0;
sk = sk3D(size(angio),Graph,name_sk,res,flag_downsampled,0);
nii = MRIread(nii_name1);
nii.vol = sk;
MRIwrite(nii,name_sk,'uchar');
%% MIP50 sk
sk = logical(sk);
ds = 50;
% sz = size(sk);
% 
% for i = 1:ds
%     I = do_downsample(sk(:,:,i:end-1-ds+i),[1,1,ds],'max');
%     if i==1; MIP = I;
%     else;    MIP = cat(2,MIP,I); 
%     end
% end
% 
% MIP = reshape(MIP,sz(1),sz(2),sz(3)-ds);
% I_last = min(sk(:,:,end-ds+1:end),[],3);
% MIP = cat(3,MIP, repmat(I_last,[1,1,ds]));

MIP = do_Z_projection(ds,sk,'max');

nii.vol = MIP;

name_sk50 = sprintf('%s/seg_%s.sk.MIP50.nii.gz',nii_dir,pad_coord);
MRIwrite(nii,name_sk50,'uchar');

%% minIP50 mus
mri = MRIread(nii_name2);

vol = single(mri.vol);
sz = size(vol);

ds = 50;
minIP = do_Z_projection(ds,vol,'min');
mri.vol = minIP;

name_I50 = sprintf('%s/I_%s.minIP50.nii.gz',nii_dir,pad_coord);
MRIwrite(mri,name_I50,'float');

%% MIP50 sk
mri = MRIread(nii_name3);

vol = single(mri.vol);
sz = size(vol);

ds = 50;
MIP = do_Z_projection(ds,vol,'max');
mri.vol = MIP;

name_I50 = sprintf('%s/sk_%s.MIP50.nii.gz',nii_dir,pad_coord);
MRIwrite(mri,name_I50,'uchar');

%% MIP50 VE
mri = MRIread(nii_name4);

vol = single(mri.vol);
sz = size(vol);

ds = 50;
MIP = do_Z_projection(ds,vol,'max');
mri.vol = MIP;

name_I50 = sprintf('%s/VE_%s.MIP50.nii.gz',nii_dir,pad_coord);
MRIwrite(mri,name_I50,'float');