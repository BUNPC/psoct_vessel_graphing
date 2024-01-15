function [refm, t] = create_mask_v4(ref, debug)
%CREATE_MASK create/apply mask for each slice
% Create boundary mask, imerode (diamond), output the mask
% INPUTS:
%   ref (double matrix): substack of psoct c-scan of tissue volume
%   debug (bool): true = show figures
% OUTPUTS:
%   refm (uint16): stack of mask for each slice
%   t (double): array of times (seconds) to run each active contour
% Mathworks References:
%   Multithresholding: Otsu, N., "A Threshold Selection Method from Gray-
%                   Level Histograms," IEEE Transactions on Systems, Man,
%                   and Cybernetics, Vol. 9, No. 1, 1979, pp. 62-66.
%   Active Contour: T. F. Chan, L. A. Vese, Active contours without edges.
%                   IEEE Transactions on Image Processing, Volume 10, Issue
%                   2, pp. 266-277, 2001.
% Overview:
%   - Use active contour to iteratively find the tissue edge for the first
%   three images each physical slice, since these have the highest
%   contrast.
%   - For each physical slice (composed of N images), project the mask from
%   the third image to the remaining images.
%   - Repeat this process for the remmaining images in the stack.

%% Initialize Parameters and Variables
% Initialize matrix for storing masked volume
refm = uint8(zeros(size(ref)));
% Create structuring element for eroding the mask
se = strel('disk',3);
% Number of iterations for active contour edge finding
nactive = 1000;
% Define array for storing timer data to track perforamnce of activecont.
t = zeros(1, size(ref, 3));
% Set the range of min/max object size (connected pixels) to retain.
min_size = 200;
max_size = size(ref,1) * size(ref,2);
range = [min_size, max_size];

%% Iterate through volume stack and find boundary mask for each slice
for ii = 1:size(ref, 3)
    %% Initialize Slice from volume + Debugging info
    s = ref(:,:,ii);
    %%% Debugging information
    fspec = 'Starting slice %i\n';
    fprintf(fspec, ii)
    
    %% Counter for b-scan within physical stack
    % Array for index within physical stack
    stackidx = mod(ii,size(ref, 3));
    % If the stack index is 0 or between 4-10, then do not compute a mask.
    % This corresponds to an image index of 4-11.
    if (1 <= stackidx) && (stackidx <= 3)
        %% Initialize mask
        %%% Multithreshold slice to find agarose pixels
        lvl = multithresh(s,13);
        % Discard the darkest portions of image (lowest threshold)
        lvl = lvl(2:end);
        % Quantize image with thresholds
        maskth = imquantize(s, lvl);
        % Set (quanta > 1) = 1 and (quanta == 1) = 0
        maskth(maskth==1) = 0;
        maskth(maskth>1) = 1;

        %%% Create histogram of voxel intensity
%         h = histogram(s);
%         v = h.Values;
%         e = h.BinEdges;
%         % Find second peak from zero intensity (this corresponds to agarose)
%         [pks,locs] = findpeaks(v);
%         pksod = sort(pks,'descend');
%         idx = pks == pksod(1);
%         p = locs(idx);
%         % Add to this to include more agarose background
%         p = p+5;
%         % Find corresponding voxel intensities from bin edges
%         agar = e(p);
%         % Remove voxels with intensity equal or less than this peak intensity
%         maskth = s;
%         maskth(maskth <= agar) = 0;


        %%% Fill holes in the mask. This is required for chan-vese method
        maskfill = imfill(maskth);
        %%% Remove islands
        mask0 = bwareafilt(logical(maskfill), range,4);
        %%% Set the method for active contour
        % The Chan-Vese method can expand or contract the mask to meet the
        % edges, whereas the "edge" method can only contract.
        ac_meth = 'Chan-vese';

        %% Iteratively find edge, erode border, remove islands
        %%% Active contour - find tissue border
        tic;
        ac = activecontour(s, mask0, nactive, ac_meth);
        t(ii) = toc;
        % Output the time for each iteration
        fspec = 'Iteration %i took %i seconds to run active contour\n';
        fprintf(fspec, ii, t(ii))
        
        %%% Erode the border + Remove islands of pixels (disjoint pixels)
        bw = imerode(ac, se);
        refm(:,:,ii) = bwareafilt(bw, range);

    else
        %%% Project mask from third slice
        % In this case, the image is between the fourth and the last image
        % in the physical slice, so the mask from the third slice will be
        % projected onto this slice.
        refm(:,:,ii) = refm(:,:,3);
    end

    %% Debugging figures
    if debug && (1 <= stackidx) && (stackidx <= 3)
        % If image 1-3 in stack, then plot these
        tmp = exist('maskth','var');
        if tmp
            figure; imshow(maskth, []); title('Multithresholded')
            figure; imshow(maskfill, []); title('Multithresholded + Filled')
        end
        % Clear maskth before next iteration
        clear maskth
        % Initial mask to feed into active contour
        figure; imshow(mask0, []);
        title('Mutlithresh, filled, islands removed',['Slice ', num2str(ii)])
        % Active contour
        figure; imshow(ac, []);
        title({'Active contour',['Slice ', num2str(ii)]})
        % Active contour + borders eroded
        figure; imshow(bw, []);
        title({'Active contour, borders eroded',['Slice ', num2str(ii)]})
        % Active contour + borders eroded + islands removed
        figure; imshow(refm(:,:,ii), []);
        title({'Active contour, border eroded, islands removed',...
            ['Slice ', num2str(ii)]})
    end

end

end

