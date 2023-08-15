% To Do
% =====
% Ithresh should be passed based on editImageThresh

function reEstimated = imView3d_CenterNodes( eventdata, centerStep, flagVisualize, flagZero, Ithresh, nLst )
global im

if eventdata==1
    ch = 2; % XY and diameter
elseif eventdata==2
    ch = 6; % straighten
end

if ~exist('flagZero')
    flagZero = [];
end
if isempty(flagZero)
    flagZero = 0;
end

if ~exist('nLst')
    nLst = [];
end
if isempty(nLst)
    nLst = 1:size(im.nodePos,1);
end

%ch = menu('Center Nodes:','XY','XY & estimate diameter','Z','Mean of surround nodes within diameter','Remove Large Angle Nodes','Towards mean of neighboring nodes','Zero diameter nodes','Cancel');

if ch==8
    reEstimated = 0;
    return;
end
reEstimated = 1;


%Ithresh = 2;


% Mean of Surround Nodes
if ch==4
    nodePos = zeros(size(im.nodePos));
    nN = size(im.nodePos,1);
    Hvox = ones(nN,1)*im.Hvox;
    for ii=1:nN
        pos = im.nodePos(ii,:);
        rhoTmp = (im.nodePos - ones(nN,1)*pos) .* Hvox;
        rho = sum(rhoTmp.*rhoTmp,2).^0.5;
        lst = find(rho<max(min(im.nodeDiam(ii)/2,20),5));
        if ~isempty(lst)
            nodePos(ii,:) = mean(im.nodePos(lst,:),1);
        else
            nodePos(ii,:) = im.nodePos(ii,:);
        end
    end
    im.nodePos = nodePos;
    return
end


% Remove Large Angle Nodes
if ch==5
    nN = size(im.nodePos,1);
    flagEdgeRemoved = 0;
    nodeFlag = ones(nN,1);
    hwait = waitbar(0,'Removing nodes with large angles');
    for ii=1:nN
        waitbar(ii/nN,hwait);
        eLst = find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii);
        if length(eLst)==2 % only operate on nodes with 2 edges
            nn = [];
            for jj=1:2  
                if im.nodeEdges(eLst(jj),1)==ii
                    pos1 = im.nodePos(im.nodeEdges(eLst(jj),1),:);
                    pos2 = im.nodePos(im.nodeEdges(eLst(jj),2),:);
                    nn(jj) = im.nodeEdges(eLst(jj),2);
                else
                    pos2 = im.nodePos(im.nodeEdges(eLst(jj),1),:);
                    pos1 = im.nodePos(im.nodeEdges(eLst(jj),2),:);
                    nn(jj) = im.nodeEdges(eLst(jj),1);
                end
                c(jj,:) = (pos2 - pos1) / norm(pos2 - pos1);
            end
            if sum(c(1,:).*c(2,:))>0.5  % remove if angle smaller than THRESH
                % only remove if connecting nodes are otherwise connected
                E = im.nodeEdges;
                E(eLst,:) = 0;
                targetFound = searchGraphForNode( E, nn(1), 5, nn(2));
                if targetFound
                    eflag = ones(size(E,1),1);
                    eflag(eLst)=0;
                    im.nodeEdges = im.nodeEdges(find(eflag==1),:);
                    flagEdgeRemoved = 1;
                    nodeFlag(ii) = 0;
                end
            end
        end
    end
    close(hwait);
    
    if flagEdgeRemoved
        im.edgeFlag = zeros(size(im.nodeEdges,1),1);

        % remove abandoned nodes
        [im.nodePos,im.nodeDiam,im.nodeDiamThetaIdx,im.nodeBC,im.nodeBCType,im.nodeType,im.nodeEdges,im.edgeFlag] = removeNodes( nodeFlag, im.nodePos, im.nodeDiam, im.nodeDiamThetaIdx, im.nodeBC, im.nodeBCType, im.nodeType, im.nodeEdges );

        reEstimated = 2;
    else
        reEstimated = 0;
    end
    return
end


% Move towards mean of neighboring nodes
if ch==6
    nN = size(im.nodePos,1);
    nodePos = im.nodePos;
    hwait = waitbar(0,'Moving nodes towards mean of neighboring nodes');
    for iii= 1:length(nLst) %1:nN
        ii = nLst(iii);
        waitbar(iii/length(nLst),hwait)
        eLst = find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii);
        nLst2 = setdiff(unique(im.nodeEdges(eLst,:)), ii);
        if length(nLst2)>1
            pos0 = max(im.nodePos(ii,:),1);
            posC = mean(im.nodePos(nLst2,:),1);
            posN = pos0 + (posC-pos0) / max(norm(posC-pos0),1);
%            if im.I(round(pos0(2)),round(pos0(1)),round(pos0(3)))>=Ithresh
                if im.I(round(posN(2)),round(posN(1)),round(posN(3)))>=Ithresh
                    nodePos(ii,:) = posN;
                end
%            else
%                nodePos(ii,:) = posN;
%            end
        end
    end
    close(hwait)
    im.nodePos = nodePos;
    reEstimated = 1;
    return
end


% Center Nodes in XY or Z
if ch<4
    minDiam = 5;
    maxDiam = 90;


    nN = size(im.nodePos,1);
    
    if ch==1 | ch==2
        pIdx = [1 2];
        nnX = im.nX;
        nnY = im.nY;
    elseif ch==3
        pIdx = [3 2];
        nnX = im.nZ;
        nnY = im.nY;
    end

    nx = ceil(maxDiam / im.Hvox(pIdx(1)));
    ny = ceil(maxDiam / im.Hvox(pIdx(2)));
    %nz = ceil(maxDiam / im.Hvox(3));
    hx = im.Hvox(pIdx(1));
    hy = im.Hvox(pIdx(2));
    %hz = im.Hvox(3);


    [nr,nc,ns] = size(im.I);
    if ch==1 | ch==2
        [nr,nc] = size(im.I(:,:,1));
    elseif ch==3
        [nr,nc] = size(squeeze(im.I(1,:,:)));
    end

    % create mapping of different angular cuts
    iTheta = 0;
    if ch==1 | ch==2
        thetaLst = [0:pi/10:pi-.1];
    elseif ch==3
        thetaLst = 0;
    end
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


    % Loop over nodes
    if ~flagVisualize
        hwait = waitbar(0,'Centering...');
    end
    for iii=1:length(nLst) %1:size(im.nodePos,1)
        ii = nLst(iii);
        if ~flagVisualize
            waitbar(iii/length(nLst),hwait);
        end

        eLst = find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii);
        
        % if centering on Z, verify that mean edge direction is not in z
        if ch==3
            eLst = find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii);
            c3 = 0;
            for jj=1:length(eLst)
                pos1 = im.nodePos(im.nodeEdges(eLst(jj),1),:);
                pos2 = im.nodePos(im.nodeEdges(eLst(jj),2),:);
                c3 = c3 + abs((pos1(3)-pos2(3)) / norm(pos1-pos2));
            end
            c3 = c3 / length(eLst);
        else
            c3 = 0;
        end
        
        % if XY centering or (z centering and c3<Thresh) then center
        if c3<0.5
            pos=round(im.nodePos(ii,:));

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


            I1 = Id(lineMap)<Ithresh;
            I2 = [];
            for ii2=1:nTheta
                I2(:,ii2) = imfill(I1(:,ii2), round(nRho/2), [0 1 0;0 1 0;0 1 0]) - I1(:,ii2);
            end

            [diam,iTheta] = min(sum(I2,1));
            % If diam = 0 the try to extrapolate out a few voxels to find the
            % vessel
            if diam == 0
                for iOff = 1:3  % this will only extrapolate +/- 3 voxels
                    for ii2=1:nTheta
                        I2(:,ii2) = imfill(I1(:,ii2), round(nRho/2)+iOff, [0 1 0;0 1 0;0 1 0]) - I1(:,ii2);
                    end
                    lstGT0 = find(sum(I2,1)>0);
                    if ~isempty(lstGT0)
                        [diam,iThetaTmp] = min(sum(I2(:,lstGT0),1));
                        iTheta = lstGT0(iThetaTmp);
                        break
                    end
                    for ii2=1:nTheta
                        I2(:,ii2) = imfill(I1(:,ii2), round(nRho/2)-iOff, [0 1 0;0 1 0;0 1 0]) - I1(:,ii2);
                    end
                    lstGT0 = find(sum(I2,1)>0);
                    if ~isempty(lstGT0)
                        [diam,iThetaTmp] = min(sum(I2(:,lstGT0),1));
                        iTheta = lstGT0(iThetaTmp);
                        break
                    end
                end
            end
            if ch==3
                iTheta = 1;
            end
            theta = thetaLst(iTheta);
            lstDiam = find(I2(:,iTheta)==1);

            if ~isempty(lstDiam)
                dR = mean(lstDiam) - maxDiam;
            else
                dR = 0;
            end
            if centerStep==1
                dR = sign(dR)*min(abs(dR),1);
            end
            %%%%%%%%%%%%%%%
            % check if this moves node off of 1:nx or 1:ny
            % if so then set dR = 0
            foo = (im.nodePos(ii,pIdx(1))+cos(theta)*dR/hx);
            if foo<1 | foo>nnX
                dR = 0;
            end
            foo = (im.nodePos(ii,pIdx(2))+sin(theta)*dR/hy);
            if foo<1 | foo>nnY
                dR = 0;
            end

            if mod(ii,10)==1 & flagVisualize
                if gcf~=12
                    figure(12);
                end
                subplot(1,3,1);
                if ch==1 | ch==2
                    imagesc(im.I(:,:,min(pos(3),ns)));
                elseif ch==3
                    imagesc(squeeze(im.I(:,pos(1),:)));
                end
                axis image
                ht=text(double(pos(pIdx(1))),double(pos(pIdx(2))),'x');
                set(ht,'color','k');
                set(ht,'fontweight','bold');
                set(ht,'horizontalalignment','center')

                Id2 = Id;
                Id2(lineMap(lstDiam,iTheta)) = 0;
                subplot(1,3,2)
                imagesc(Id2,[0 32])
                axis image
                ht=text(nx+1,ny+1,'x');
                set(ht,'color','k');
                set(ht,'fontweight','bold');
                set(ht,'horizontalalignment','center')
                ht=text(nx+1+cos(theta)*dR/hx,ny+1+sin(theta)*dR/hy,'x');
                set(ht,'color','m');
                set(ht,'fontweight','bold');
                set(ht,'horizontalalignment','center')
                title( sprintf('Node %d',ii) )

                cm = jet(32);
                cm(1,:) = [1 1 1];
                colormap(cm)

                subplot(1,3,3)
                imagesc(I2+(1-I1))

                drawnow
            end

            if ch==2
                im.nodeDiamEst(ii) = diam;
                im.nodeDiam(ii) = diam;
                im.nodeDiamThetaIdx(ii) = iTheta;
            end
            im.nodePos(ii,pIdx) = im.nodePos(ii,pIdx) + [cos(theta)/hx sin(theta)/hy]*dR;
        end % end check on c3 thresh for z centering
    end % end loop over nodes
    if ~flagVisualize
        close(hwait);
    end

end % end handling ch<4

if ~isfield(im,'nodeDiam')
    return;
end

lst = find(im.nodeDiam==0);

if (length(lst)<size(im.nodeDiam,2) | ch==2) & ~isempty(lst)
    if ~flagZero
        ch = 1;
%        ch = menu( sprintf('%d nodes with diam=0. Center them between neighboring nodes?',length(lst)),'Yes','No');
    else
        ch=1;
    end
    
    if ch==2
        return
    end
    
    for ii=1:length(lst)
        iN = lst(ii);
        lstE = find(im.nodeEdges(:,1)==iN | im.nodeEdges(:,2)==iN);
        lstN = setdiff( unique( im.nodeEdges(lstE,:) ), iN );
        
        if lstN>=2
            im.nodePos(iN,:) = mean( im.nodePos(lstN,:), 1 );
        end
    end
end

        


