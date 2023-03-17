img = TIFF2MAT('ves_seg_grayscale.tif');
img_new = single(zeros(size(img)));

for i=1:size(img,3)
    temp = img(:,:,i);
    p = prctile(temp(:),99.5);
    temp(temp(:)<=p) = 0;
    temp(temp(:)>p) = 255;
    img_new(:,:,i) = temp;
end

% remove short unconnected components
CC = bwconncomp(imbinarize(img_new));
I_seg = img_new;
for uuu = 1:length(CC.PixelIdxList)
    if length(CC.PixelIdxList{uuu}) < 30    % 30 for default
        I_seg(CC.PixelIdxList{uuu}) = 0;
    end
end
MAT2TIFF(uint8(I_seg),'ves_seg.tif');

