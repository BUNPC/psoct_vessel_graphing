function [volm, t] = create_mask_v1(vol, debug, NSLOTS)
%CREATE_MASK create/apply mask for each slice
% Create boundary mask, imerode (diamond), output the mask
% INPUTS:
%   vol (double matrix): psoct c-scan of tissue volume
% OUTPUTS:
%   volm (uint16): stack of mask for each slice
%   t (double): array of times (seconds) to run each active contour
% Mathworks References:
%   Multithresholding: Otsu, N., "A Threshold Selection Method from Gray-
%                   Level Histograms," IEEE Transactions on Systems, Man,
%                   and Cybernetics, Vol. 9, No. 1, 1979, pp. 62-66.
%   Active Contour: T. F. Chan, L. A. Vese, Active contours without edges.
%                   IEEE Transactions on Image Processing, Volume 10, Issue
%                   2, pp. 266-277, 2001.

%% Initialize Parameters and Variables
% Initialize matrix for storing masked volume
volm = uint16(zeros(size(vol)));
% Create structuring element for eroding the mask
se = strel('disk',15);
% Number of iterations for active contour edge finding
nactive = 5e3;
% Define array for storing timer data to track perforamnce of activecont.
t = zeros(1, size(vol, 3));

%% Iterate through volume stack and find boundary mask for each slice
parfor (ii = 1:size(vol, 3), NSLOTS)
    %% Initialize Slice
    % Extract slice from volume and segmentation
    s = vol(:,:,ii);
    % Set the range of min/max object size (connected pixels) to retain.
    min_size = 9e4;
    max_size = size(s,1) * size(s,2);
    range = [min_size, max_size];
      
    %% Create starting mask for active contour
    %%% Create histogram of voxel intensity
    h = histogram(s);
    v = h.Values;
    e = h.BinEdges;
    % Find second peak from zero intensity (this corresponds to agarose)
    [pks,locs] = findpeaks(v);
    pksod = sort(pks,'descend');
    idx = pks == pksod(1);
    p = locs(idx);
    % Add to this to include more agarose background
    p = p+10;
    % Find corresponding voxel intensities from bin edges
    agar = e(p);
    % Remove voxels with intensity equal or less than this peak intensity
    mask0 = s;
    mask0(mask0 <= agar) = 0;
    % Fill mask
    maskfill = imfill(mask0);
    % Remove islands
    mask = bwareafilt(logical(maskfill), range,4);

    %%% Multithreshold slice to find agarose pixels
%     lvl = multithresh(s,20);
%     % Discard the darkest portions of image (lowest threshold)
%     lvl = lvl(2:end);
%     % Quantize image with thresholds
%     maskth = imquantize(s, lvl);
%     % Set (quanta > 1) = 1 and (quanta == 1) = 0
%     maskth(maskth==1) = 0;
%     maskth(maskth>1) = 1;
%     % Fill mask
%     maskfill = imfill(maskth);
%     % Remove islands
%     mask = bwareafilt(logical(maskfill), range,4);
    
    %% Active contour - find tissue border
    tic;
    ac = activecontour(s, mask, nactive);
    t(ii) = toc;
    % Output the time for each iteration
    fspec = 'Iteration %i took %i seconds to run active contour\n';
    fprintf(fspec, ii, t(ii))
    
    %% Erode the border
    bw = imerode(ac, se);
    
    %% Remove islands of pixels
    % Set the range of min/max object size (connected pixels) to retain.
    volm(:,:,ii) = bwareafilt(bw, range);

    %% Debugging figures
    if debug
        tmp = exist('mask0','var');
        if tmp
            figure; imshow(mask0, []); title('Histogram Third Peak Removed')
        else
            figure; imshow(maskth, []); title('Multithresholded')
        end
        figure; imshow(maskfill, []); title('Multithresholded + Filled')
        figure; imshow(mask, []); title('Multithresh, fill, rm islands')
        % Active contour
        figure; imshow(ac, []);
        title({'Active contour',[num2str(nactive),' Iterations']})
        % Active contour + borders eroded
        figure; imshow(bw, []);
        title({'Active contour + border eroded',[num2str(nactive),' Iterations']})
        % Active contour + borders eroded + islands removed
        figure; imshow(volm(:,:,ii), []);
        title({'Active contour + border eroded + islands removed',[num2str(nactive),' Iterations']})
    end

end

end

