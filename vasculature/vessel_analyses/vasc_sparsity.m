function [vsp] = vasc_sparsity(angio, mask, vox_dim)
%VASC_SPARSITY compute the distance of each tissue voxel to a blood vessel
%   For the tissue voxel under consideration, search for vessels within a
%   search radius. If there are no vessels within the search radius, then
%   increase the radius until it contains one. The linear indices of the
%   tissue voxels in this radius are stored in the array "vidcs"
%
%   Then, find the Euclidean distance between each vessel voxel and the
%   tissue voxel under consideration. This is stored in the array "d".
%
%   Finally, find the minimum distance of all distances, and store this in
%   the vessel sparsity matrix vsp(ii).
%
%   INPUTS:
%       seg (logical matrix): segmentation matrix
%       mask (logical matrix): tissue mask for region
%       vox_dim (double array): dimensions of 
%   OUTPUTS:
%       vs (double array): Euclidean distance of each non-vessel voxel from
%           the nearest blood vessel (microns)
%

%% Initialize paralell pool (if not yet running)
% Set # threads = # cores for job
nslots = str2num(getenv('NSLOTS'));
maxNumCompThreads(nslots);

% Check to see if we already have a parpool, if not create one with
% our desired parameters
myPool = gcp('nocreate');
if isempty(myPool)
    myCluster=parcluster('local');
    % Ensure multiple parpool jobs don't overwrite other's temp files
    myCluster.JobStorageLocation = getenv('TMPDIR');
    myPool = parpool(myCluster, nslots);
end

%% Calculate vascular sparsity

%%% Identify indices of non-vessel voxels in the masked OCT matrix
% Remove the vessel indices from tissue indices
tiss_idx = mask;
tiss_idx(angio) = [];
% Create array of non-vessel tissue indices (non-zero)
tiss_idx = find(tiss_idx);
% Convert tissue linear indices to matrix subscripts
[y,x,z] = ind2sub(size(angio), tiss_idx);

% Initialize vessel sparsity array
vsp = zeros(length(tiss_idx),1);

% Create local copy of voxel dimensions
vox_y = vox_dim(1);
vox_x = vox_dim(2);
vox_z = vox_dim(3);

%%% Specify the index ranges for each parfor worker
% Divide number of indices by number of workers
n_idcs = floor(length(tiss_idx) / nslots);
% Set the index range to the number of indices per worker
opts = parforOptions(myPool,...
                    'RangePartitionMethod','fixed',...
                    'SubrangeSize',n_idcs);

% Iterate over tissue indices
parfor (ii = 1:length(tiss_idx), opts)
    %%% Expand search radius (ROI) until it contains a vessel voxel
    % vessel voxel indices within search radius
    vidcs = [];
    % Search radius length (voxels)
    r = 50;
    while isempty(vidcs)
        % Increment the search radius
        r = r + 10;
        % Define ROI for finding vessel voxels
        yi = max(1, y(ii) - r);
        yf = min(y(ii) + r, size(angio,1));
        xi = max(1, x(ii) - r);
        xf = min(x(ii) + r, size(angio,2));
        zi = max(1, z(ii) - r);
        zf = min(z(ii) + r, size(angio,3));
        % Find vessels in ROI
        vidcs = find(angio(yi:yf,xi:xf,zi:zf));
    end
    
    %% Find euclidean distance between ROI center and all vessel voxels
    
    %%% Convert matrix subscripts to offset in microns
    % Convert linear indices of vessel voxels to matrix subscripts
    [vy,vx,vz] = ind2sub(size(angio(yi:yf,xi:xf,zi:zf)),vidcs);
    % Scale vessel voxel matrix subscripts by voxel dimensions (microns)
    vy = vy .* vox_y;
    vx = vx .* vox_x;
    vz = vz .* vox_z;
    
    %%% Calculate offset of the voxel in the center of the ROI
    y0 = ((size(angio(yi:yf,xi:xf,zi:zf),1) - 1)/2 + 1) .* vox_y;
    x0 = ((size(angio(yi:yf,xi:xf,zi:zf),2) - 1)/2 + 1) .* vox_x;
    z0 = ((size(angio(yi:yf,xi:xf,zi:zf),3) - 1)/2 + 1) .* vox_z;
    
    %%% Euclidean dist. b/w center voxel and vessel voxels in ROi
    % Initialize distance matrix
    d = zeros(length(vidcs),1);
    % Iterate over each vessel voxel
    for j = 1:length(d)
        % Compute distance between vessel and tissue
        d(j) = sqrt((vy(j)-y0).^2 + (vx(j)-x0).^2 + (vz(j)-z0).^2);
    end
    % Find minimum of all Euclidean distances and store
    vsp(ii) = min(d(:));
end


end