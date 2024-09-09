function ld = length_density(data, t_mask)
%% Calculate the length density in 1/microns^2
% Find all non-zero voxels in tissue volume, convert to logical, take sum
% to identify total number of tissue voxels, then convert to cubic microns.
% Then, divide the total length of all vessels by the total volume.
%
% INPUTS:
%   data (struct): data structure of graph
%   t_mask (uint8): tissue mask
% OUTPUTS:
%   ld (double): length density

% Retrieve voxel dimensions (x,y,z) in microns
vox = data.Graph.vox;
% Calculate volume of a single voxel (cubic microns)
vox_vol = vox(1) .* vox(2) .* vox(3);

% Calculate volume of tissue from non-zero voxels in tissue mask (cubic
% microns)
t_vol = sum(logical(t_mask(:))) .* vox_vol;

% Calculate total length of vasculature (microns)
len = sum(data.Graph.segInfo.segLen_um);

% Calculate length density from sum of all segments and total tissue vol
ld = len ./ t_vol;

% Convert length density from um/um^3 -> mm/mm^3
ld = ld .* 1e6;

end

    