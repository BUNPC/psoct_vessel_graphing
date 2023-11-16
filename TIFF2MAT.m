function dat=TIFF2MAT(tiffname)
%% Convert tiff image into matrix
% author: Jiarui Yang
% 02/21/20

% return tiff structure, one element per image
tiff_info = imfinfo(tiffname);
% read in first image in stack
dat = imread(tiffname, 1) ;
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread(tiffname, ii);
    dat = cat(3 , dat, temp_tiff);
end
end