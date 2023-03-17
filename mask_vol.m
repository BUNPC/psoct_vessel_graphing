function mask=mask_vol(vol)
% function that creates a binray mask for tissue region
% Author: Jiarui Yang
% 02/27/20
mask=false(size(vol));
for z=1:size(vol,3)
    Img=squeeze(vol(:,:,z));
    % Gaussian filtering
    sigma=3;     % scale parameter in Gaussian kernel
    G=fspecial('gaussian',15,sigma);
    Img_smooth=conv2(Img,G,'same');
    Img_smooth=uint16(Img_smooth);
    % find level to binraize image
    level=graythresh(Img_smooth);
    BW=imbinarize(Img_smooth,level);
    BW_filled = imfill(BW,'holes');
    if sum(BW_filled(:))<0.9*size(BW_filled,1)*size(BW_filled,2)
        boundaries = bwboundaries(BW_filled);
        len=zeros(1,length(boundaries));
        for i=1:length(boundaries)
            b=boundaries{i};
            len(i)=length(b);
        end
        [~, loc]=max(len);
        cord=boundaries{loc};
        for p=1:size(cord,1)
            mask(cord(p,1),cord(p,2),z)=true;
        end
        mask(:,:,z)=imfill(squeeze(mask(:,:,z)),'holes');
    else
        mask(:,:,z)=mask(:,:,z-1);
    end
end

end