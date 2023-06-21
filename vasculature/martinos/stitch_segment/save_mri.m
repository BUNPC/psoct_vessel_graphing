function save_mri(I, name, res, datatype)
% save_mri save the volume
% INPUTS:
%       I (matrix): volume
%       name (string): output filename
%       res (array): [x,y,z] dimensions
%       datatype (string): output datatype

% Make Nifti and header
colres = res(2); 
rowres = res(1); 
sliceres = res(3); 
mri.volres = [res(1) res(2) res(3)];
mri.xsize = rowres;
mri.ysize = colres;
mri.zsize = sliceres;
a = diag([-colres rowres sliceres 1]);
mri.vox2ras0 = a;
mri.volsize = size(I);
mri.vol = I;

% Save output
disp('Started saving volume.')
MRIwrite(mri,name,datatype);
disp('Finished saving volume.');

end