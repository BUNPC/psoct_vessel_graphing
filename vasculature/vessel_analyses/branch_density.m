function bd = branch_density(data, t_mask)
%% Calculate the branch density in 1/microns^3
% Find all non-zero voxels in tissue volume, convert to logical, take sum
% to identify total number of tissue voxels, then convert to cubic microns.
% Then, divide the total number of branchesby the total volume.
%
% INPUTS:
%   data (struct): data structure of graph
%   t_mask (uint8): tissue mask
% OUTPUTS:
%   bd (double): branch density

% Retrieve voxel dimensions (x,y,z)
vox = data.Graph.vox;
% Calculate volume of a single voxel
vox_vol = vox(1) .* vox(2) .* vox(3);
% Calculate volume of tissue from the non-zero voxels in tissue mask
t_vol = sum(logical(t_mask(:))) .* vox_vol;
% Retrieve number of branch points
nb = data.Graph.nB;
% Calculate branch density
bd = sum(nb) / t_vol;
% Convert branch density to cubic millimeters
bd = bd .* 1e9;
end