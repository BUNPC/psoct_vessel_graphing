function I = loadImage( im )

if strcmp(im.filetype,'img')
    filenm = sprintf( '%s.img',im.filenm );
    fid = fopen( filenm, 'rb' );
    n = fread( fid, 3, 'int');
    nx = n(1); ny = n(2); nz = n(3);
    I = fread( fid, nx*ny*nz, 'float' );
    I = reshape(I, [ny nx nz] );
    fclose(fid);

elseif strcmp(im.filetype,'tif') | strcmp(im.filetype,'tiff')
    finfo = imfinfo(sprintf('%s.%s',im.filenm,im.filetype));
    ny = finfo(1).Height;
    nx = finfo(1).Width;
    nz = length(finfo);
    I = zeros(ny,nx,nz);
    for ii=1:nz
        I(:,:,ii) = imread(sprintf('%s.%s',im.filenm,im.filetype),ii);
    end
    I = I / max(I(:));
        
elseif strcmp(im.filetype,'mat')
    load(im.filenm,'I');
    I = I / max(I(:));
    
end

