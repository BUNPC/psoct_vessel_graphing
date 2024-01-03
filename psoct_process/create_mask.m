function [volm, segm] = create_mask(vol, seg, ref)
%CREATE_MASK create/apply mask for each slice
% Create boundary mask, imerode (diamond), output the mask
% INPUTS:
%   vol (double matrix): psoct b-scan of tissue volume
%   seg (uint8): segmentation volume
% OUTPUTS:
%   volm (uint16): masked volume
%   segm (uint8): masked segmentation

%% Initialize
% Initialize matrix for storing masked volume
volm = uint16(zeros(size(vol)));
segm = uint16(zeros(size(seg)));
% Create structuring element for eroding the mask
se = strel('diamond',40);

%% Iterate through volume stack and find boundary mask for each slice
for ii = 1:size(vol, 3)
    % Extract slice from volume and segmentation
    slice = vol(:,:,ii);
    
    %%% Boundary Mask
    % Create boundary mask for current slice
    bmask = boundarymask(slice);
    % Erode the mask
    bmask = imerode(bmask,se);
    % Fill holes of mask. This is to avoid false negative where the lumen
    % is classified as background and removed.
    bmask = imfill(bmask, "holes");
    % Range of object size to keep
    range = [1e4, 1e8];
    bmask = bwareafilt(bmask, range);

    %%% Pia mask
    % Calculate differences between pixel intensities
    w = graydiffweight(slice, ref, 'GrayDifferenceCutoff', 10000);
    % Threshold the difference map
    wth = log(w);
    wth(wth<4) = 0;    
    % Dilate the segmentation
    se = strel('disk',5);
    wth = imdilate(wth, se);    
    % Keep components within range of connections
    wth = logical(wth);
    range = [1e4, 1e8];
    bw = bwareafilt(wth, range);
    % Dilate again to expand mask
    bw = imdilate(bw, se);
    % Invert mask
    piamask = ~bw;
    % Again, only keep components within range of connections
    piamask = bwareafilt(piamask, range);

    %%% Combine background mask and pia mask --> make joint mask
    mask = piamask .* bmask;
    % Dilate mask
    se = strel('square',5);
    mask = imerode(mask, se);

    %%% Apply mask to volume and segmentation
    % Apply mask to tissue slice
    slice_me = slice .* uint16(mask);
    % Apply mask to segmentation slice
    seg_me = uint8(mask) .* uint8(seg(:,:,ii));
    % Add the masked slice to the new stack
    volm(:,:,ii) = slice_me;
    segm(:,:,ii) = seg_me;

    %%% Debugging plots
%     figure; imshow(piamask);
%     figure; imshow(mask);
%     figure; imshow(slice_me);

end

end

