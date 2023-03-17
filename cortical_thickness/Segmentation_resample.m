addpath /autofs/cluster/octdata2/users/Hui/ProcessHK001_cerebellum_20200618
basepath = '/autofs/cluster/octdata2/users/Hui/ProcessHK001_cerebellum_20200618/';
segv = load_nifti([basepath 'Stack_nii/seg.nii']);
data = segv.vol;
mol = single(data==85);
gra = single(data==170);
pix_x = segv.pixdim(2);
pix_z = segv.pixdim(3);
pix_y = segv.pixdim(4);
downsample = 3;
[Xq, Yq, Zq] = meshgrid(linspace(1,size(data,2),size(data,2)*round(pix_z/pix_x/downsample)), linspace(1, size(data,1), round(size(data,1)/3)), linspace(1, size(data,3), round(size(data,3)/3)));

mol_up = interp3(mol, Xq, Yq, Zq);
gra_up = interp3(gra, Xq, Yq, Zq);


resolution = pix_x*downsample * ones(1,3);

% mkdir('Thickness');
name_mol = [basepath 'Thickness/seg_resample_mol.nii'];
MakeNii(name_mol, mol_up, resolution);

name_gra = [basepath 'Thickness/seg_resample_gra.nii'];
MakeNii(name_gra, gra_up, resolution);



