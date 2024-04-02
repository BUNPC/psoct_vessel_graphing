function overlay_vol_seg(vol, seg, color, fout, invert_color)
%overlay_vol_seg Overlay the volume with the segmentation output
% INPUTS:
%   vol (mat): tissue volume (grayscale)
%   seg (mat): segmentations of vasculature (binary)
%   color (string): color of segmentation in overlaid image
%   fout (string): full filepath for output .TIF

%% Convert volume back to non-inverted (black background)
if invert_color
    vol = imcomplement(vol);
end

%% Overlay segmentation and probaiblity
nvol = size(vol,3);
for j=1:nvol
    tmp = imoverlay(vol(:,:,j), seg(:,:,j), color);
    if j==1
        ov = tmp;
    else
        ov = cat(4,ov,tmp);
    end
end

rgb_stack_2_tif(ov, fout)

