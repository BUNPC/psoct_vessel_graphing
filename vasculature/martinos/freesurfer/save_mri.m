function save_mri(I, name, res, datatype, permuteflag)
% SAVE_MRI convert matrix to NII
% INPUTS:
%   I (matrix): input volume
%   name (string): output filename. options include:
%           1. MGH file. Eg, f.mgh or f.mgz
%           2. BHDR file Eg, f.bhdr. Result will be written to a bfloat
%                   volume, eg f_000.bfloat.
%           3. NIFIT file, Eg, f.nii, f.nii.gz (uncompressed and compressed)
%   res (double array, 1x3): resolution of each dimension [x,y,z]
%   datatype (string): uchar, short, int, float, double, ushort, uint
%              Only applies to nifti. Setting datatype to '' implies float.
%           
%   permuteflag (0/1): 0 = retain order of [x,y,z,t]
%                      1 = swap x,y, ie [y,x,z,t]
    
%% Make the header file for the MRI
disp(' - making hdr...');
rowres = res(1);
colres = res(2);
sliceres = res(3);
mri.volres = [res(1) res(2) res(3)];
mri.xsize = rowres;
mri.ysize = colres;
mri.zsize = sliceres;
a = diag([-colres rowres sliceres 1]);
mri.vox2ras0 = a;
mri.volsize = size(I);
mri.vol = I;

%% Write the MRI
disp('- writing MRI -')
MRIwrite(mri,name,datatype,permuteflag);
disp(' - done - ');

end