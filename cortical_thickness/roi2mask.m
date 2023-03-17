function roi2mask(input_name, roi_name, output_name)
    addpath('E:\Jiarui\Data\code\');
    addpath('E:\Jiarui\Data\code\cortical_thickness');
    img = TIFF2MAT(input_name);
    mask = zeros(size(img));
    roi = ReadImageJROI(roi_name);
    roi_region = ROIs2Regions(roi,[size(img,1), size(img,2)]);
    for i = 1:length(roi_region.PixelIdxList)
        slice = zeros(size(img,2), size(img,1));
        pxl_idx_list = roi_region.PixelIdxList{i};
        slice(pxl_idx_list) = 255;
        mask(:,:,i) = transpose(slice);
    end
    MAT2TIFF(single(mask),output_name);
end