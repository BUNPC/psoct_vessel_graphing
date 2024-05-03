% Creates the line-segment masks and stores them in a sparse format.
%
% L = makeLineSegments(radius, angleStep, upSampleScale)
%
% radius:        Radius of the sphere containing the line segments. It
%                should be a few voxels larger than the expected maximum
%                thickness.
% angleStep:     The directions of the line segments will be 'angleStep'
%                degrees apart. (Default: 10)
% upSampleScale: Subvoxel resolution. Increasing this value softens the
%                edges of the masks through partial voluming. (Default: 5)
% useParpool:    If true, the parallel pool will be used. (Default: false.
%                If a parallel pool exists, it will automatically be used.)
%
% L:             Created line-segment masks.
%
% Reference:
% I. Aganj, G. Sapiro, N. Parikshak, S. K. Madsen, and P. Thompson,
% "Measurement of cortical thickness from MRI by minimum line integrals on
% soft-classified tissue," Human Brain Mapping, vol. 30, no. 10,
% pp. 3188-3199, 2009. http://doi.org/10.1002/hbm.20740
%
% See also:   MLI_thickness, EXAMPLE.m

% Codes by Iman Aganj.

function L = makeLineSegments(radius, angleStep, upSampleScale, useParpool)

if ~exist('angleStep', 'var')
    angleStep = 10;
end
if ~exist('upSampleScale', 'var')
    upSampleScale = 5;
end
if ~isempty(gcp('nocreate')) || (exist('useParpool', 'var') && useParpool)
    nCores = getfield(gcp, 'NumWorkers');
else
    nCores = 0;
end

for theta = (0:angleStep:89) * pi/180
    if theta == 0
        angles = [0; 0];
    else
        phi = (0:(angleStep/sin(theta)):359) * pi/180;
        angles = [angles, [theta*ones(1,length(phi)); phi]];
    end
end
direc = [bsxfun(@times, sin(angles(1,:)), [cos(angles(2,:)); sin(angles(2,:))]); cos(angles(1,:))];
nDirecs = size(direc,2);
[X,Y,Z] = ndgrid(-(radius+2):(radius+2));
R2 = X.^2 + Y.^2 + Z.^2;
[XI,YI,ZI] = meshgrid(-(radius+2)*upSampleScale-(upSampleScale-1)/2:(radius+2)*upSampleScale+(upSampleScale-1)/2);
XI = XI/upSampleScale; YI = YI/upSampleScale; ZI = ZI/upSampleScale;
RI2 = XI.^2 + YI.^2 + ZI.^2;
L = cell(1,nDirecs);
parfor (n = 1:nDirecs, nCores)
    IP = direc(1,n)*XI + direc(2,n)*YI + direc(3,n)*ZI;
    C = (IP>=0) & (IP<=radius) & (RI2 - IP.^2 <= .25);
    L0 = zeros(2*radius+5, 2*radius+5, 2*radius+5);
    for i = 1:2*radius+5
        for j = 1:2*radius+5
            for k = 1:2*radius+5
                CP = C(upSampleScale*(i-1)+1:upSampleScale*i,upSampleScale*(j-1)+1:upSampleScale*j,upSampleScale*(k-1)+1:upSampleScale*k);
                L0(i,j,k) = sum(CP(:));
            end
        end
    end
    L0 = L0 * radius / sum(L0(:));
    M = L0(:)>0;
    [~, ind] = sort(R2(M));
    Lt = [X(M) Y(M) Z(M) L0(M)];
    L{n} = Lt(ind, :);
end
