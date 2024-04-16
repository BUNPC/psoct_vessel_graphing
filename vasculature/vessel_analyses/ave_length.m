function [um_len, vox_len] = ave_length(data)
%% Calculates mean/average length of vessels in um and voxels 
% Input: graph .mat file 
    um_len = mean(data.Graph.segInfo.segLen_um);
    vox_len = mean(data.Graph.segInfo.segLen);
end