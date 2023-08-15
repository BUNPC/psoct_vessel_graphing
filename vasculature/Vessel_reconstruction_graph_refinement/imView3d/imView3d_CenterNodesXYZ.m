% To Do
% Image thresh from GUI
% use imfill to mask out unconnected vessels
%   but what if node not connected to any vessels
% option to pass nB

function imView3d_CenterNodesXYZ( centerStep1vox, flagVisualize, IThresh, nLst )
global im


if ~exist('nLst')
    nLst = [];
end
if isempty(nLst)
    nLst = 1:size(im.nodePos,1);
end

I = im.I;
%IThresh = 2;

nNodes = size(im.nodePos,1);
[ny,nx,nz] = size(I);
nodePos = im.nodePos;


if ~flagVisualize
    hwait = waitbar(0,'Centering XYZ...');
end

nB = zeros(nNodes,1);
for ii=1:nNodes
    nB(ii) = length(find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii));
end

for jNode=1:length(nLst) %1:nNodes
    iNode = nLst(jNode);

    if ~flagVisualize
        waitbar(jNode/length(nLst),hwait);
    end

    if nB(iNode)<=2
        [lstE,lstC] = find(im.nodeEdges==iNode);

        pos1 = im.nodePos(iNode,:);
        pos2 = im.nodePos(im.nodeEdges(lstE(1),mod(lstC(1),2)+1),:);
        if length(lstE)==2
            pos0 = im.nodePos(im.nodeEdges(lstE(2),mod(lstC(2),2)+1),:);
        else
            pos0 = [];
        end

        r = norm(pos2-pos1);
        if r>0
            theta = acos((pos2(3)-pos1(3))/r);
            rho = norm(pos2(1:2)-pos1(1:2));
            phi = acos((pos2(1)-pos1(1))/rho);

            if~isempty(pos0)
                r = norm(pos1-pos0);
                if r>0
                    theta2 = acos((pos1(3)-pos0(3))/r);
                    rho = norm(pos1(1:2)-pos0(1:2));
                    phi2 = acos((pos1(1)-pos0(1))/rho);

                    while (phi-phi2)>3.14159
                        phi2 = phi2 + 2*3.14159;
                    end
                    theta = (theta + theta2)/2;
                    phi = (phi+phi2)/2;
                end
            end

            xLst = [-10:10];
            yLst = [-10:10];
            Isub = zeros(length(yLst),length(xLst));
            for ix = 1:length(xLst)
                for iy = 1:length(yLst)
                    dpos = [xLst(ix) yLst(iy) 0]';
                    dpos = [[-cos(theta) 0 sin(theta)];[0 1 0];[-sin(theta) 0 cos(theta)]] * dpos;
                    dpos = [[cos(phi) -sin(phi) 0];[sin(phi) cos(phi) 0];[0 0 1]] * dpos;
                    iix = max(min(round(pos1(1)+dpos(1)),nx),1);
                    iiy = max(min(round(pos1(2)+dpos(2)),ny),1);
                    iiz = max(min(round(pos1(3)+dpos(3)),nz),1);
                    Isub(iy,ix) = I( iiy, iix, iiz );
                end
            end

            Isub = Isub .* (Isub >= IThresh); % I could add a continuity condition here
            % using imfill
            IsubSum = max(sum(Isub(:)),1)+eps;
            [xx,yy] = meshgrid(xLst,yLst);
            posM(1) = sum(xx(:).*Isub(:))/IsubSum;
            posM(2) = sum(yy(:).*Isub(:))/IsubSum;

            if flagVisualize
                figure(1)
                subplot(2,1,1)
                imagesc(I(:,:,round(pos1(3))))
                xlim([-20 20]+pos1(1));
                ylim([-20 20]+pos1(2));
                subplot(2,1,2)
                imagesc(xLst,yLst,Isub)
                ht = text(0,0,'o');
                ht = text(posM(1),posM(2),'x');
            end


            if centerStep1vox
                posM = posM / (max(norm(posM),1)+eps);
            end

            dpos = [posM(1) posM(2) 0]';
            dpos = [[-cos(theta) 0 sin(theta)];[0 1 0];[-sin(theta) 0 cos(theta)]] * dpos;
            dpos = [[cos(phi) -sin(phi) 0];[sin(phi) cos(phi) 0];[0 0 1]] * dpos;
            iix = max(min(pos1(1)+dpos(1),nx),1);
            iiy = max(min(pos1(2)+dpos(2),ny),1);
            iiz = max(min(pos1(3)+dpos(3),nz),1);
            nodePos(iNode,:) = [iix iiy iiz];
            %        pause(0.1)
        end % r>0
    end  % end if nB(iNode)<=2
end % end loop on iNode

im.nodePos = nodePos;

close(hwait)
