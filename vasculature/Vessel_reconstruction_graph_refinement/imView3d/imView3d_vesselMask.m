function imView3d_vesselMask( )
global im

nx = im.nX;
ny = im.nY;
nz = im.nZ;

if ~isfield(im,'Hvox')
    im.Hvox = 1;
    hwait = waitbar(0,'Please enter the voxel size!');
    while length(im.Hvox~=3)
        im.Hvox = str2num(input('Enter voxel size [x y z]?','s'));
        if length(im.Hvox~=3)
            disp('Please enter 3 lengths!');
        end
    end
    close(hwait);
end

if isfield(im,'edgeVel')
    edgeVel = abs(im.edgeVel);
    maxEdgeVel = max(edgeVel);
    maxEdgeVel = 10e3;
    edgeFlow = abs(im.edgeFlow);
    maxEdgeFlow = max(edgeFlow);
    nodePressure = im.nodePressure;
    maxPressure = max(im.nodePressure(:));
else
    edgeVel = zeros(size(im.nodeEdges,1),1);
    maxEdgeVel = 1;
    edgeFlow = zeros(size(im.nodeEdges,1),1);
    maxEdgeFlow = 1;
    nodePressure = zeros(size(im.nodePos,1),1);
    maxPressure = 1;
end

updateFlag = 0;
if ~isfield(im,'nodeTypeUpdated')
    im.nodeTypeUpdated = ones(length(im.nodeType),1);
elseif ~isempty(find(im.nodeTypeUpdated==1))
    ch = menu('Updated mask for all nodes or just modified node?','All','Modified');
    if ch==1
        im.nodeTypeUpdated = ones(length(im.nodeType),1);
    else
        updateFlag = 1;
    end
else
    im.nodeTypeUpdated = ones(length(im.nodeType),1);
end


nEdges = size(im.nodeEdges,1);
nNodes = size(im.nodePos,1);

hwait = waitbar(0,'Allocating memory for mask');
Vm = uint8(zeros(ny,nx,nz));
waitbar(.25,hwait);
if isfield(im,'edgeVel')
    Vvel = uint8(zeros(ny,nx,nz));
    waitbar(.5,hwait);
    Vflow = uint8(zeros(ny,nx,nz));
    waitbar(.75,hwait);
    Vpres = uint8(zeros(ny,nx,nz));
end
close(hwait)

if updateFlag==1
    Vm = im.Vm;
    if isfield(im,'edgeVel')
        Vvel = im.Vvel;
        Vflow = im.Vflow;
        Vpres = im.Vpres;
    end
end

%hWait = waitbar( 0, sprintf('Masking vessels...\nMake sure diameters\nand velocities\nhave been estimated.'));
hWait = waitbar( 0, sprintf('Masking vessels...'));
nodeMasked = zeros(nNodes,1);
for iE=1:nEdges
    waitbar(iE/nEdges,hWait);
    i1 = im.nodeEdges(iE,1);
    i2 = im.nodeEdges(iE,2);

    if (im.nodeTypeUpdated(i1) | im.nodeTypeUpdated(i2)) & (nodeMasked(i1)==0 | nodeMasked(i2)==0)
        r1 = im.nodeDiam(i1)/2;
        if r1<1
            r1=2;
        elseif r1>10
            r1=10;
        end
        r2 = im.nodeDiam(i2)/2;
        if r2<1
            r2=2;
        elseif r2>10
            r2=10;
        end
        p1 = round(im.nodePos(i1,:));
        p2 = round(im.nodePos(i2,:));
        d12 = norm(p2-p1);
        dxyz = (p2-p1)/d12;
        rd = (r2-r1)/d12;

        if im.nodeType(i1)>0
            nType = im.nodeType(i1)+1;
        else
            nType = 1;
        end
        
        lst = find(sum((im.nodePos-ones(nNodes,1)*p1).^2,2).^0.5<(r1) );
        nodeMasked(lst) = 1;
        lst = find(sum((im.nodePos-ones(nNodes,1)*p2).^2,2).^0.5<(r2) );
        nodeMasked(lst) = 1;

        p = p1;
        r = r1;
        flag = 1;
        while norm(round(p)-p2)~=0 || flag
            if norm(round(p)-p2)==0
                flag = 0;
            end
            pr = round(p);
            rTmp = min(r,10);
            rx = ceil(rTmp/im.Hvox(1));
            ry = ceil(rTmp/im.Hvox(2));
            rz = ceil(rTmp/im.Hvox(3));
            for iX = -rx:+rx
                for iY = -ry:+ry
                    for iZ = -rz:+rz
                        if norm([iX*im.Hvox(1) iY*im.Hvox(2) iZ*im.Hvox(3)])<=r
                            iix = min(max(pr(1)+iX,1),nx);
                            iiy = min(max(pr(2)+iY,1),ny);
                            iiz = min(max(pr(3)+iZ,1),nz);
                            Vm(iiy,iix,iiz) = 1 + nType;
                            if isfield(im,'edgeVel')
                                Vvel(iiy,iix,iiz) = round(32*min(edgeVel(iE)/maxEdgeVel,1));
                                Vflow(iiy,iix,iiz) = round(32*min(edgeFlow(iE)/maxEdgeFlow,1));
                                Vpres(iiy,iix,iiz) = round(32*min(nodePressure(i1)/maxPressure,1));
                            end
                        end
                    end
                end
            end

            %         for iX = max(pr(1)-rx,1):min(pr(1)+rx,nx)
            %             for iY = max(pr(2)-ry,1):min(pr(2)+ry,ny)
            %                 for iZ = max(pr(3)-rz,1):min(pr(3)+rz,nz)
            %                     if norm([iX-pr(1) iY-pr(2) iZ-pr(3)])<=r
            %                         Vm(iY,iX,iZ) = 1 + nType;
            %                     end
            %                 end
            %             end
            %         end
            if flag
                p = p + dxyz;
                r = r + rd;
            end
        end
    end % end of check if node type updated
end
close(hWait);

im.nodeTypeUpdated = zeros(length(im.nodeType),1);

im.Vm = uint8(Vm);
if isfield(im,'edgeVel')
    im.Vvel = uint8(Vvel);
    im.Vflow = uint8(Vflow);
    im.Vpres = uint8(Vpres);
end