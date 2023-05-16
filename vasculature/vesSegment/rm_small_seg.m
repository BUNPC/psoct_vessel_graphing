function [I_seg, cnt, seg_total] = rm_small_seg(prob_mat, min_conn)
%rm_small_seg Remove small segments of vessels
%   Identify segments of connected voxels with fewer than "min_conn"
%   connected voxels. Remove these small segments because they are likely
%   false positives.
%
%   INPUTS:
%       prob_mat (matrix): output of the frangi filter
%       min_conn (double): minimum number of connections. Segments with
%                       fewer will be removed.
%   OUTPUTS:
%       I_seg (matrix): probability matrix after removing small connections.
%       cnt (double): number of removed segments
%       cc.NumObjects (double): total number of segments prior to removal

% Run built-in matlab function to perform connectivity analysis
cc = bwconncomp(prob_mat);
seg_total = cc.NumObjects;
% Define output variable
I_seg = prob_mat;
% Counter to track number of vessels that are removed
cnt = 0;
for ii = 1:length(cc.PixelIdxList)
    if length(cc.PixelIdxList{ii}) < min_conn
        I_seg(cc.PixelIdxList{ii}) = 0;
        cnt = cnt + 1;
    end
end

end