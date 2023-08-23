% To Do
% =====
% Ithresh should be passed based on editImageThresh
%
% after estimating diameter I should consider including more z-slices to
%    make sure that diameter does not change

function reEstimated = imView3d_GetDiam()
global im


ch = menu('Estimate Node Diameters.','All','New','Cancel');

if ch==3
    reEstimated = 0;
    return;
end
reEstimated = 1;


minDiam = 5;
maxDiam = 90;

Ithresh = 2;


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


for ii=1:size(im.nodePos,1)

    if ch==1 | (ch==2 & im.nodeDiam(ii)==1)
        pos=round(im.nodePos(ii,:));

        if (pos(1)-nx)>=1
            dx1 = 1;
            sx1 = pos(1)-nx;
        else
            dx1 = nx-pos(1)+2;
            sx1 = 1;
        end
        if (pos(1)+nx)<=nc
            dx2 = 2*nx+1;
            sx2 = pos(1)+nx;
        else
            dx2 = 2*nx+1-((pos(1)+nx)-nc);
            sx2 = nc;
        end

        if (pos(2)-ny)>=1
            dy1 = 1;
            sy1 = pos(2)-ny;
        else
            dy1 = ny-pos(2)+2;
            sy1 = 1;
        end
        if (pos(2)+ny)<=nr
            dy2 = 2*ny+1;
            sy2 = pos(2)+ny;
        else
            dy2 = 2*ny+1-((pos(2)+ny)-nr);
            sy2 = nr;
        end

        Id = 32*ones(2*ny+1, 2*nx+1);
        Id(dy1:dy2,dx1:dx2) = im.I(sy1:sy2,sx1:sx2,pos(3));
        Id(find(Id<Ithresh))=0;


        I1 = Id(lineMap)<Ithresh;
        I2 = [];
        for ii2=1:nTheta
            I2(:,ii2) = imfill(I1(:,ii2), round(nRho/2), [0 1 0;0 1 0;0 1 0]) - I1(:,ii2);
        end

        [diam,iTheta] = min(sum(I2,1));
        theta = thetaLst(iTheta);
        lstDiam = find(I2(:,iTheta)==1);
        
        if ~isempty(lstDiam)
            dR = mean(lstDiam) - maxDiam;
        else
            dR = 0;
        end
        %%%%%%%%%%%%%%%
        % check if this moves node off of 1:nx or 1:ny
        % if so then set dR = 0
        foo = (im.nodePos(ii,1)+cos(theta)*dR/hx);
        if foo<1 | foo>im.nX
            dR = 0;
        end
        foo = (im.nodePos(ii,2)+sin(theta)*dR/hy);
        if foo<1 | foo>im.nY
            dR = 0;
        end

        figure(12);
        subplot(1,3,1);
        imagesc(im.I(:,:,pos(3)));
        axis image
        ht=text(double(pos(1)),double(pos(2)),'x');
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
%         ht=text(nx+1+cos(theta)*dR/hx,ny+1+sin(theta)*dR/hy,'x');
%         set(ht,'color','m');
%         set(ht,'fontweight','bold');
%         set(ht,'horizontalalignment','center')
        title( sprintf('Node %d,  Diam = %d um',ii, diam) )

        cm = jet(32);
        cm(1,:) = [1 1 1];
        colormap(cm)

        subplot(1,3,3)
        imagesc(I2+(1-I1))

        drawnow

        im.nodeDiamEst(ii) = diam;
        im.nodeDiam(ii) = diam;
%        im.nodePos(ii,:) = im.nodePos(ii,:) + [cos(theta)/hx sin(theta)/hy 0]*dR;
    end
end



