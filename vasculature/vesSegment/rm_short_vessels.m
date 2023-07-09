function [seg_out] = rm_short_vessels(seg_in, voxmin)
%% Remove small disconnected segments via connectivity analysis
% Remove segments that are composed of fewer than 30 voxels.
%
% INPUTS:
%       seg_in (matrix): output of Frangi filter
%       voxmin (int): minimum number of voxels connected to classify vessel
%
% OUTPUTS:
%       seg_out (matrix): parsed output of Frangi filter
%       cnt (int): number of vessels removed during connectivity analysis

% Run built-in matlab function to perform connectivity analysis
cc = bwconncomp(seg_in);
% Define output variable
seg_out = seg_in;
% Counter to track number of vessels that are removed
cnt = 0;
for ii = 1:length(cc.PixelIdxList)
    if length(cc.PixelIdxList{ii}) < voxmin
        seg_out(cc.PixelIdxList{ii}) = 0;
        cnt = cnt + 1;
    end
end

sprintf('Segments before connectivity analysis = %d', cc.NumObjects)
sprintf('Segments removed during connectivity analysis = %d', cnt)
sprintf('Segments after connectivity analysis = %d', (cc.NumObjects - cnt))

end

