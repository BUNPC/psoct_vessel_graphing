% Computes the tissue thickness from a probability image, using minimum
% line integrals.
%
% MLI_thickness_MEX.c must be compiled using "mex MLI_thickness_MEX.c".
% Alternatively, the compiled files can be downloaded from:
% www.nitrc.org/projects/thickness
% 
% [Thickness, SkeletonDistance] = MLI_thickness(I, L, param)
%  
% I:       3D image representing a soft segmentation (probability map) of
%          the tissue, with values between 0 and 1.
% L:       Line segments, created by makeLineSegments.m
% param:   Optional parameter structure containing the following fields:
%       param.threshStop:  The threshold for the first stopping criterion
%          (see the paper). During line integration, if the image value has
%          been below this threshold for param.numVoxStop consecutive
%          voxels, the integration will stop. (Default: 0.3)
%       param.numVoxStop:  See param.threshStop above. (Default: 1)
%       param.numVoxValey: The number of voxels for the second stopping
%          criterion (see the paper). During line integration, if the image
%          has been decreasing for this number of successive voxels and
%          then increasing for an additional such number of voxels, the
%          integration will stop. (Default: 0; i.e., disabled).
%       param.useParpool:  If 'true', a parallel pool will be used for
%          speedup. (Default: false)
% 
% Thickness:         Thickness image.
% SkeletonDistance:  The distance-transform of the skeleton. E.g.,
%                    (SkeletonDistance<1) is a 2-voxel-wide skeleton mask.
% 
% Reference:
% I. Aganj, G. Sapiro, N. Parikshak, S. K. Madsen, and P. Thompson,
% "Measurement of cortical thickness from MRI by minimum line integrals on
% soft-classified tissue," Human Brain Mapping, vol. 30, no. 10,
% pp. 3188-3199, 2009. http://doi.org/10.1002/hbm.20740
%
% See also:   makeLineSegments, EXAMPLE.m

% Codes by Iman Aganj.

function [Thickness, SkeletonDistance] = MLI_thickness(I, L, param)

if ~exist('param', 'var')
    param = [];
end
if ~isfield(param, 'threshStop')
    param.threshStop = .3;
end
if ~isfield(param, 'numVoxStop')
    param.numVoxStop = 1;
end
if ~isfield(param, 'numVoxValey')
    param.numVoxValey = 0;  % Disables the second stopping criterion.
end
if ~isfield(param, 'useParpool')
    param.useParpool = false;
end

try
    if param.useParpool
        nCores = getfield(gcp, 'NumWorkers');
        Thickness = zeros(size(I));
        oneSidedThickness = Thickness;
        parParts = round(linspace(1, numel(I)+1, nCores+1));
        parfor n = 1:nCores
            [Thickness_partial, ~, oneSidedThickness_partial] = MLI_thickness_MEX(double(I), L, double(param.numVoxValey), double(param.threshStop), double(param.numVoxStop), [parParts(n), parParts(n+1)-1]);
            Thickness = Thickness + Thickness_partial;
            oneSidedThickness = oneSidedThickness + oneSidedThickness_partial;
        end
    else
        [Thickness, ~, oneSidedThickness] = MLI_thickness_MEX(double(I), L, double(param.numVoxValey), double(param.threshStop), double(param.numVoxStop)); % The ignored (second) output parameter is the index of the optimally chosen line direction.
    end
catch ME
    if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
        warning(sprintf('The compiled mex file was not found. Use "mex MLI_thickness_MEX.c" to compile MLI_thickness_MEX.c.\nAlternatively, you can download the compiled files from: www.nitrc.org/projects/thickness'))
    end
    rethrow(ME)
end

if nargout > 1
    SkeletonDistance = abs(.5*Thickness - oneSidedThickness);
end
