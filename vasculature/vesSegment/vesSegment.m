function [w, I_seg] = vesSegment(I, gsigma, gsize, pmin, voxmin)
% This function performs 3D vessel enhancemnt, filtering, and thesholding
% to segment the vessels from the surrounding tissue.
%
% INPUTS:
%   I (matrix) - inverted 3D angiogram (white = vessel, non-white = tissue).
%   gsigma (array): vector of standard deviation values for gaussian filter
%   gsize (vector): Size of the 3D gaussian kernel (voxels). Must be
%                       either a single positive, odd integer, or a
%                       3-element array of positive, odd integers. If a
%                       single positive, odd integer is specified (Q), then
%                       the gsize will be [Q, Q, Q].
%   pmin - threshold to determine which voxel belongs to a vessel. This is
%           applied to the probability matrix from the output of the frangi
%           filter. The threshold belongs to [0, 1].
%   voxmin - minimum number of connected voxels to classify as vessels
%
% OUTPUTS:
%   w (matrix) - probability of voxel belonging to vessel
%   I_seg (matrix) - binary matrix. 1 = vessel. 0 = non-vessel tissue. This
%               matrix is the result of thresholding the matrix w.
%{
Copyright (c) 2011, Zhang Jiang
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
%}
    
%% Initialization
% Frangi filter parameters
alpha = 0.25;    
gamma12 = 0.5;
gamma23 = 0.5;

% Initialize eigenvalue matrix and output volume
[k,l,m] = size(I);
Lambda123 = zeros(k,l,m);
w = zeros(k,l,m);

%% Vessel enhancement filtering (frangi filter)
h = waitbar(0,'Performing vessel enhancement');
for ii = 1:length(gsigma)
    waitbar((ii-1)/length(gsigma));
    
    % Apply Hessian3D filter (Gaussian kernel) then calculate 2nd order
    % gradients (approximation of 2nd order derivatives of image).
    [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(I, gsigma(ii), gsize(ii));
    
    %%% Calculate eigenvalues (L1, L2, L3) of image Hessians
    % Call executable depending on operating system
    if ispc
        [L1, L2, L3, ~, ~, ~] =...
            eig3volume_windows(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
    else
        [L1, L2, L3, ~, ~, ~] =...
            eig3volume_linux(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
    end
    
    SortL = sort([L1(:)'; L2(:)'; L3(:)'],'descend');
    L1 = reshape(SortL(1,:),size(L1));
    L2 = reshape(SortL(2,:),size(L2));
    L3 = reshape(SortL(3,:),size(L3));
    
    %%% Find eigenvalues that match conditions for frangi vessel
    % NOTE: the conditions for Lambda123 can be reduced when gamma12 = gamma23 = 0.5;
    %       refer to page 4 of the 1998 Sato et al. paper
    %       "Three-dimensional multi-scale line filter for segmentation..."
    % First condition
    idx = find((L3 < L2) & (L2 < L1) & (L1 <= 0));
    if (gamma12 == gamma23) && (gamma12 == 0.5)
        Lambda123(idx) = sqrt( abs(L3(idx)) .* (abs(L2(idx) + L1(idx))));
    else
        Lambda123(idx) = abs(L3(idx)).*(L2(idx)./L3(idx)).^gamma23.*...
                        (1+L1(idx)./abs(L2(idx))).^gamma12;
    end
    % Second condition
    idx = find((L3 < L2) & (L2 < 0) & (0 < L1) & (L1 < abs(L2)./alpha));
    if (gamma12 == gamma23) && (gamma12 == 0.5)
        Lambda123(idx) = sqrt( abs(L3(idx)) .* (abs(L2(idx) - alpha.*L1(idx))));
    else
        Lambda123(idx) = abs(L3(idx)).*(L2(idx)./L3(idx)).^gamma23.*...
                        (1-alpha*L1(idx)./abs(L2(idx))).^gamma12;
    end
    %%% Take maximum likelihood values for each vesselness measure (w)
    % If there is only one sigma value, then set w equal to the first
    % output of Lamba123
    if ii == 1
        w = Lambda123;
    else
        w = max(w, Lambda123);
    end
end
close(h);

%%% normalize the data (probably normalizing slice wise might be good?)
w = (w-min(w(:))) / (max(w(:))-min(w(:)));

%% Apply threshold to probability matrix.
% The matrix w represents the likelihood each voxel belongs to a vessel.
% The threshold (pmin) determines the cutoff for this likelihood.
% (w < pmin) = 0 (non-vessel)
% (w < pmin) = 1 (vessel)
w(w<pmin) = 0;
w(w>=pmin) = 1;

%% Remove small disconnected segments via connectivity analysis
% Remove segments that are composed of fewer than 30 voxels.
I_seg = rm_short_vessels(w_thresh, voxmin);
disp('Done segmentation')



