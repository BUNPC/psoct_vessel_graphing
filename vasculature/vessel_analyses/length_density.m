function ld = length_density(data, t_mask)
%% Calculate the length density in 1/microns^2
% Find all non-zero voxels in tissue volume, convert to logical, take sum
% to identify total number of tissue voxels, then convert to voxels. Then,
% divide the total length of all vessels by the total volume.
%
% INPUTS:
%   data (struct): data structure of graph
%   t_mask (uint8): tissue mask
% OUTPUTS:
%   ld (double): length density

    vox = data.Graph.vox;
    vox_vol = vox(1) .* vox(2) .* vox(3);
    t_vol = sum(logical(t_mask(:))) .* vox_vol;
    ld = sum(data.Graph.segInfo.segLen_um) ./ t_vol;
end

    