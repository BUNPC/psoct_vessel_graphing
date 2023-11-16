function imView3d_erodeCenterLine( Ithresh )
global im

% Erosion center line
%Ithresh = 2;
I = im.I >= Ithresh;

%%
% create mapping of different angular cuts
minDiam = 5;
maxDiam = 90;


nN = size(im.nodePos,1);

pIdx = [1 2];
nnX = im.nX;
nnY = im.nY;

nx = ceil(maxDiam / im.Hvox(pIdx(1)));
ny = ceil(maxDiam / im.Hvox(pIdx(2)));
%nz = ceil(maxDiam / im.Hvox(3));
hx = im.Hvox(pIdx(1));
hy = im.Hvox(pIdx(2));
%hz = im.Hvox(3);

hx = 1; % i assume that hx=hy and so might as well go by voxel steo rather than 1 um step
hy = 1;
nx = maxDiam;
ny = maxDiam;

[nr,nc,ns] = size(im.I);
[nr,nc] = size(im.I(:,:,1));

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





% first run a center xy to make sure we dealt with all nodes at I<Ithresh
% any nodes at I<Ithresh will remain there
[nIx,nIy,nIz] = size(I);
nN = size(im.nodePos,1);
lstIdx = min(round(im.nodePos(:,2)),nIy) + (min(round(im.nodePos(:,1)),nIx)-1)*nIy + (min(round(im.nodePos(:,3)),nIz)-1)*nIx*nIy;
nLst = find( I( lstIdx )==0 );
nodeFixed = zeros(nN,1);
nodeFixed(nLst) = 1;

% LOOP
I2 = I;
iter = 1;
nI2 = nr*nc*ns;
nI2gt0 = nI2;
nMoved = 0;
while length(find(nodeFixed==1))<nN & nI2gt0>0 
    hwait = waitbar(0,sprintf('Erosion centerline iteration %d... %d fixed out of %d... %d gt0 of %d... Moved %d',iter,length(find(nodeFixed==1)),nN,nI2gt0,nI2,nMoved));

    I = I2;
    % erode image in XY planes
    for ii=1:ns
        I2(:,:,ii) = imerode(I(:,:,ii),strel('square',3));
    end
    nI2gt0 = length(find(I2>0));

    % find nodes that are at I<Ithresh
    lstIdx = min(round(im.nodePos(:,2)),nIy) + (min(round(im.nodePos(:,1)),nIx)-1)*nIy + (min(round(im.nodePos(:,3)),nIz)-1)*nIx*nIy;
%    lstIdx = round(im.nodePos(:,2)) + (round(im.nodePos(:,1))-1)*nIy + (round(im.nodePos(:,3))-1)*nIx*nIy;
    nLst = find( I2( lstIdx )==0 );
    nLst = setdiff( nLst, find(nodeFixed==1) );
%    disp(sprintf('Checking %d nodes + fixed %d = %d', length(nLst),length(find(nodeFixed==1)),length(find(nodeFixed==1))+length(nLst)) )
    
    %%
    % center nLst
    nMoved = 0;
    for ii=1:length(nLst)
        waitbar(ii/length(nLst),hwait);
        iN = nLst(ii);

        pos=round(im.nodePos(iN,:));

        % extract the sub image 
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
        
        % check z slices up and down 
        for iZ = [0 -1 1 -2 2]
            Id = 32*ones(2*ny+1, 2*nx+1);
            Id(dy1:dy2,dx1:dx2) = I2(sy1:sy2,sx1:sx2,min(max(pos(3)+iZ,1),ns));
            %    Id(find(Id<Ithresh))=0;


            iTheta = im.nodeDiamThetaIdx(iN);
            Id1 = Id(lineMap(:,iTheta)) == 0;%<Ithresh;
            Id2 = [];
            Id2 = imfill(Id1, round(nRho/2), [0 1 0;0 1 0;0 1 0]) - Id1;

            [diam] = min(sum(Id2,1));
            % If diam = 0 then try to extrapolate out a few voxels to find the
            % vessel
            iOff = 0;
            if diam == 0
                for iOff = 1:2  % this will only extrapolate +/- 2 voxels
                    Id2 = imfill(Id1, round(nRho/2)+iOff, [0 1 0;0 1 0;0 1 0]) - Id1;
                    lstGT0 = find(sum(Id2,1)>0);
                    if ~isempty(lstGT0)
                        [diam] = min(sum(Id2(:,lstGT0),1));
                        break
                    end
                    Id2 = imfill(Id1, round(nRho/2)-iOff, [0 1 0;0 1 0;0 1 0]) - Id1;
                    lstGT0 = find(sum(Id2,1)>0);
                    if ~isempty(lstGT0)
                        [diam] = min(sum(Id2(:,lstGT0),1));
                        break
                    end
                end
            end
            if diam ~= 0
                im.nodePos(iN,3) = im.nodePos(iN,3) + iZ;
                break
            end
        end
        

        % update node XY position
        theta = thetaLst(iTheta);
        lstDiam = find(Id2==1);

        if ~isempty(lstDiam)
            dR = mean(lstDiam) - maxDiam - 1;
        else
            dR = 0;
        end
        dR = sign(dR) * iOff;
        foo = (im.nodePos(iN,pIdx(1))+cos(theta)*dR/hx);
        if foo<1 | foo>nnX
            dR = 0;
        end
        foo = (im.nodePos(iN,pIdx(2))+sin(theta)*dR/hy);
        if foo<1 | foo>nnY
            dR = 0;
        end
        if dR~=0
            nMoved = nMoved + 1;
            im.nodePos(iN,pIdx) = im.nodePos(iN,pIdx) + [cos(theta)/hx sin(theta)/hy]*dR;
        end
        
    end

    % those nodes that are still at I<Ithresh are now fixed
    lstIdx = min(round(im.nodePos(:,2)),nIy) + (min(round(im.nodePos(:,1)),nIx)-1)*nIy + (min(round(im.nodePos(:,3)),nIz)-1)*nIx*nIy;
%    lstIdx = round(im.nodePos(:,2)) + (round(im.nodePos(:,1))-1)*nIy + (round(im.nodePos(:,3))-1)*nIx*nIy;
    nLst = find( I2( lstIdx )==0 );
    nodeFixed(nLst) = 1;


    % loop until all nodes fixed
    iter = iter + 1;
    close(hwait)
end
















    