function nodeDiam = calc_diameter(angio, nodes, Ithresh, vox_dim)
% Calculate diameter for each node position. 
%   For each node, create a z-axis projection of the segmentation. Then,
%   calculate the diameter at each node.
%
%   INPUTS:
%       angio (logical matrix) - PSOCT volume or segmentation volume
%       nodes (double matrix) - node position from graph (in voxels)
%       edges (double matrix) - edge indices from graph
%       Ithresh () - approximate threshold of the segmentation
%       vox_dim () - voxel size (microns)
%          
%   OUTPUTS:
%       nodeDiam - Diameter at each node (microns). If the voxels are
%       isotropic in X-Y then this will be a 1D array.
%       Otherwise, it will be a 2D matrix where each column contains the
%       minimum and maximum diameter for each node.
%   
%   This is modified GetDiam function from imView3d
%   Modified by Sreekanth Kura - skura@bu.edu
%   Further modified by Mack Hyman - mhyman@bu

%% Initialize variables
% Set minimum and maximum diameters
minDiam = 2;
maxDiam = 30;

% initialize vector for storing diameter at each node
nodeDiam = zeros(1,size(nodes,1));

% Calculate 
nx = ceil(maxDiam / vox_dim(1));
ny = ceil(maxDiam / vox_dim(2));

% Voxel dimensions in X,Y
hx = vox_dim(1);
hy = vox_dim(2);

% number of voxels in X,Y,Z dimensions
[nr,nc,ns] = size(angio);

%% Create mapping of different angular cuts
thetaLst = 0:pi/10:pi-.1;
nTheta = length(thetaLst);
rhoLst = -maxDiam:1:maxDiam;
nRho = length(rhoLst);
lineMap = zeros(nRho,nTheta);
for iTheta = 1:nTheta
    for iRho = 1:nRho
        xx = rhoLst(iRho) * cos(thetaLst(iTheta));
        yy = rhoLst(iRho) * sin(thetaLst(iTheta));

        ix = round((xx+maxDiam)/hx);
        iy = round((yy+maxDiam)/hy);

        lineMap(iRho,iTheta) = (2*ny+1)*ix + iy + 1;
    end
end

%% Calculate diameter for each node
nNodes = size(nodes,1);
for ii = 1:nNodes    
    % Retrieve position of nodes (voxel graph index)
    pos = round(nodes(ii,:));
    
    %%% Set the range for X,Y in the Id and angio
    % X Dimension - lower bound
    if (pos(1)-nx)>=1
        dx1 = 1;
        sx1 = pos(1)-nx;
    else
        dx1 = nx-pos(1)+2;
        sx1 = 1;
    end
    % X Dimension - upper bound
    if (pos(1)+nx)<=nc
        dx2 = 2*nx+1;
        sx2 = pos(1)+nx;
    else
        dx2 = 2*nx+1-((pos(1)+nx)-nc);
        sx2 = nc;
    end
    
    % Y Dimension - lower bound
    if (pos(2)-ny)>=1
        dy1 = 1;
        sy1 = pos(2)-ny;
    else
        dy1 = ny-pos(2)+2;
        sy1 = 1;
    end
    % Y Dimension - upper bound
    if (pos(2)+ny)<=nr
        dy2 = 2*ny+1;
        sy2 = pos(2)+ny;
    else
        dy2 = 2*ny+1-((pos(2)+ny)-nr);
        sy2 = nr;
    end
    
    Id = 32*ones(2*ny+1, 2*nx+1);
    zrange = min(max(pos(3)-1,1),ns):min(max(pos(3)+1,1),ns);
    Id(dy1:dy2,dx1:dx2) = max(angio(sy1:sy2,sx1:sx2,zrange),[],3);
    
    % Set Id below threshold equal to zero
    Id(Id<Ithresh) = 0;
    

    I1 = Id(lineMap)<Ithresh;
    I2 = zeros(size(I1));
    for ii2=1:nTheta
        I2(:,ii2) = imfill(I1(:,ii2), round(nRho/2),...
                            [0 1 0;0 1 0;0 1 0]) - I1(:,ii2);
    end
    
    [diam,~] = min(sum(I2,1));
    % If diam = 0, extrapolate out to find vessel
    if diam == 0
        % Extrapolate +/- 3 voxels
        for iOff = 1:3
            for ii2=1:nTheta
                I2(:,ii2) = imfill(I1(:,ii2), round(nRho/2)+iOff,...
                                    [0 1 0;0 1 0;0 1 0]) - I1(:,ii2);
            end
            lstGT0 = find(sum(I2,1)>0);
            if ~isempty(lstGT0)
                [diam,~] = min(sum(I2(:,lstGT0),1));
                break
            end
            
            for ii2=1:nTheta
                I2(:,ii2) = imfill(I1(:,ii2), round(nRho/2)-iOff,...
                                    [0 1 0;0 1 0;0 1 0]) - I1(:,ii2);
            end
            lstGT0 = find(sum(I2,1)>0);
            if ~isempty(lstGT0)
                [diam,~] = min(sum(I2(:,lstGT0),1));
                break
            end
        end
    end
    
    %make sure we respect threshold
    if diam<minDiam
        diam=minDiam;
    elseif diam>maxDiam
        diam=maxDiam;
    end
    
    % Set the outputs
    nodeDiam(ii) = diam;
end

%%% Scale diameter by the voxel dimension
% If isotropic, then just use x-dimension
if hx == hy
    nodeDiam = nodeDiam .* hx;
% otherwise, output a range
else
    nodeDiam_min = nodeDiam .* min([hx,hy]);
    nodeDiam_max = nodeDiam .* max([hx,hy]);
    nodeDiam = [nodeDiam_min; nodeDiam_max];
end