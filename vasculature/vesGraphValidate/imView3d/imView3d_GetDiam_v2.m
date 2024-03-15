% To Do
% =====
% Ithresh should be passed based on editImageThresh
%
% after estimating diameter I should consider including more z-slices to
%    make sure that diameter does not change
%
% This currently only takes a single 

function reEstimated = imView3d_GetDiam_v2( Ithresh )
global im


if isfield(im,'estDiamFlag')  % flag for running diamter code automatically
    if im.estDiamFlag==1
        ch = 1;
    else
        ch = menu('Estimate Node Diameters.','All','New','Cancel');
    end
else
    ch = menu('Estimate Node Diameters.','All','New','Cancel');
end

if ch==3
    reEstimated = 0;
    return;
end
reEstimated = 1;


minDiam = 0.5;
maxDiam = 50;

%Ithresh = 2;


nx = ceil(maxDiam / im.Hvox(1));
ny = ceil(maxDiam / im.Hvox(2));
nz = ceil(maxDiam / im.Hvox(3));
hx = im.Hvox(1);
hy = im.Hvox(2);
hz = im.Hvox(3);


[nr,nc,ns] = size(im.I);


% create mapping of different angular cuts
iTheta = 0;
thetaLst = [0:pi/10:pi-.1];
nTheta = length(thetaLst);
rhoLst = [-maxDiam:1:maxDiam];
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

pIdx = [1 2];

hwait = waitbar(0,'Get Diameter...');
nNodes = size(im.nodePos,1);
for ii=1:size(im.nodePos,1)
    
    if ch==1 | (ch==2 & im.nodeDiam(ii)==1)
        
        waitbar(ii/nNodes,hwait);
        
        eLst = find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii);
        
        pos=round(im.nodePos(ii,:));
        pos(pos == 0) = 1; %- Added by Allen Alfadhel on 30-April-2020
        
        %% What does this section do?
        if (pos(pIdx(1))-nx)>=1
            dx1 = 1;
            sx1 = pos(pIdx(1))-nx;
        else
            dx1 = nx-pos(pIdx(1))+2;
            sx1 = 1;
        end

        if (pos(pIdx(1))+nx)<=nc
            dx2 = 2*nx+1;
            sx2 = pos(pIdx(1))+nx;
        else
            dx2 = 2*nx+1-((pos(pIdx(1))+nx)-nc);
            sx2 = nc;
        end
        
        if (pos(pIdx(2))-ny)>=1
            dy1 = 1;
            sy1 = pos(pIdx(2))-ny;
        else
            dy1 = ny-pos(pIdx(2))+2;
            sy1 = 1;
        end
        
        if (pos(pIdx(2))+ny)<=nr
            dy2 = 2*ny+1;
            sy2 = pos(pIdx(2))+ny;
        else
            dy2 = 2*ny+1-((pos(pIdx(2))+ny)-nr);
            sy2 = nr;
        end
        
        Id = 32*ones(2*ny+1, 2*nx+1);
        if ch==1 | ch==2
            Id(dy1:dy2,dx1:dx2) = im.I(sy1:sy2,sx1:sx2,min(pos(3),ns));
        elseif ch==3
            Id(dy1:dy2,dx1:dx2) = squeeze(im.I(sy1:sy2,pos(1),sx1:sx2));
        end
        Id(find(Id<Ithresh))=0;
        
        %% What does this section do?
        I1 = Id(lineMap) < Ithresh;
        I2 = [];
        for j=1:nTheta
            I2(:,j) = imfill(I1(:,j), round(nRho/2), [0 1 0;0 1 0;0 1 0]) - I1(:,j);
        end
        
        [diam,iTheta] = min(sum(I2,1));
        % If diam = 0 the try to extrapolate out a few voxels to find the
        % vessel
        if diam == 0
            for iOff = 1:3  % this will only extrapolate +/- 3 voxels
                for j=1:nTheta
                    I2(:,j) = imfill(I1(:,j), round(nRho/2)+iOff, [0 1 0;0 1 0;0 1 0]) - I1(:,j);
                end
                lstGT0 = find(sum(I2,1)>0);
                if ~isempty(lstGT0)
                    [diam,iThetaTmp] = min(sum(I2(:,lstGT0),1));
                    iTheta = lstGT0(iThetaTmp);
                    break
                end
                for j=1:nTheta
                    I2(:,j) = imfill(I1(:,j), round(nRho/2)-iOff, [0 1 0;0 1 0;0 1 0]) - I1(:,j);
                end
                lstGT0 = find(sum(I2,1)>0);
                if ~isempty(lstGT0)
                    [diam,iThetaTmp] = min(sum(I2(:,lstGT0),1));
                    iTheta = lstGT0(iThetaTmp);
                    break
                end
            end
        end
        
        %%% If the measured diameter exceeds limits then set to limits
        if diam<minDiam
            diam=minDiam;
        elseif diam>maxDiam
            diam=maxDiam;
        end
        
        im.nodeDiamEst(ii) = diam;
        im.nodeDiam(ii) = diam;
        im.nodeDiamThetaIdx(ii) = iTheta;        
        
    end
end

if length(im.nodeDiamEst)>nNodes
    im.nodeDiamEst(nNodes+1:end) = [];
end

%%% Calculate median diameter for each segment
for sg = 1:length(im.segDiam) %- added by Allen 12-May-20
    im.segDiam(sg) = median(im.nodeDiam(im.nodeSegN == sg));
end


close(hwait)