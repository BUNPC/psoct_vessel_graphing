function [I_VE,I_seg] = vesSegment(I, sigma, thres)
% This function performs 3D vessel enhancemnt, filtering, and thesholding
% to segment the vessels from the surrounding tissue.
%
% INPUTS:
%   I - inverted 3D angiogram (white = vessel, non-white = tissue).
%   sigma - vector of standard deviation values of gaussian filter to
%           calcualte hessian matrix at each voxel
%   thres - threshold to determine which voxel belongs to a vessel. This is
%           applied to the probability matrix from the output of the frangi
%           filter.
% OUTPUTS:
%   I_VE () - likelihood of voxel belonging to vessel
%   I_seg () - binary matrix. 1 = vessel. 0 = non-vessel tissue. This
%               matrix is the result of applying the threshold to I_VE.
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
    
% set parameters
alpha = 0.25;    
gamma12 = 0.5;
gamma23 = 0.5;

% get volume size and cerate the empty voulme for output
[k,l,m] = size(I);
I_VE = zeros(k,l,m);
T_temp = zeros(k,l,m);

% vessel enhancement filtering
h = waitbar(0,'Please wait... performing vessel enhancement');
for i = 1:length(sigma)
    waitbar((i-1)/length(sigma));

    [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(I,sigma(i));
    if ispc 
        [Lambda1, Lambda2, Lambda3, ~, ~, ~] =...
            eig3volume_windows(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
    else
        [Lambda1, Lambda2, Lambda3, ~, ~, ~] =...
            eig3volume_linux(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
    end
    
    SortL = sort([Lambda1(:)'; Lambda2(:)'; Lambda3(:)'],'descend');
    Lambda1 = reshape(SortL(1,:),size(Lambda1));
    Lambda2 = reshape(SortL(2,:),size(Lambda2));
    Lambda3 = reshape(SortL(3,:),size(Lambda3));
    
    idx = find(Lambda3 < 0 & Lambda2 < 0 & Lambda1 < 0);
    T_temp(idx) = abs(Lambda3(idx)).*(Lambda2(idx)./Lambda3(idx)).^gamma23.*(1+Lambda1(idx)./abs(Lambda2(idx))).^gamma12;
    
    idx = find(Lambda3 < 0 & Lambda2 < 0 & Lambda1 > 0 & Lambda1 < abs(Lambda2)/alpha);
    T_temp(idx) = abs(Lambda3(idx)).*(Lambda2(idx)./Lambda3(idx)).^gamma23.*(1-alpha*Lambda1(idx)./abs(Lambda2(idx))).^gamma12;
    if i == 1
        I_VE = T_temp;
    else
        I_VE = max(I_VE,T_temp);
    end
end
close(h);

% normalize the data (probably normalizing slice wise might be good?)
I_VE = (I_VE-min(I_VE(:)))/(max(I_VE(:))-min(I_VE(:)));

%% Apply threshold to probability matrix.
% The matrix I_VE represents the likelihood each voxel belongs to a vessel.
% The threshold (thres) determines the cutoff for this likelihood.
% I_VE < thres = 0 (non-vessel)
% I_VE < thres = 0 (vessel vessel)
T_thres = I_VE;
T_thres(I_VE<thres) = 0;
T_thres(I_VE>=thres) = 1;

%% remove small disconnected segments via connectivity analysis
CC = bwconncomp(T_thres);
I_seg = T_thres;
for uuu = 1:length(CC.PixelIdxList)
    if length(CC.PixelIdxList{uuu}) < 30    % 30 for default
        I_seg(CC.PixelIdxList{uuu}) = 0;
    end
end