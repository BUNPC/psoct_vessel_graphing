function [im] = center_nodes_xyz(im, centerStep1vox, visualize_flag, im_thresh, nLst)
%% center_nodes_XYZ move nodes towards center line of segment
%{

INPUTS:
    im (struct): Matlab structure containing
                    angio (matrix): volumetric segmentation
                    nodes: [y,x,z] coordinates of nodes
                    edges: start/end indices of segments
    centerStep1vox ():
    visualize_flag (bool): whether or not to display updated data
    im_thresh ():
    nLst ():

OUTPUTS:
    im (struct): updated with new nodes

To Do (past comments from David):
- use imfill to mask out unconnected vessels
    - but what if node not connected to any vessels
- option to pass nB
%}


%% Initialize nLst, angio, nNodes, nodes, nB, [ny,nx,nz]

% Initialize nLst if it does not exist
if ~exist('nLst', 'var')
    nLst = [];
end
if isempty(nLst)
    nLst = 1:size(im.nodes,1);
end

% Initialize segmentation variable
angio = im.angio;

% Initialize nodes
nNodes = size(im.nodes,1);
[ny,nx,nz] = size(angio);
nodes = im.nodes;

% Initialize number of branches (nB)
nB = zeros(nNodes,1);
for ii=1:nNodes
    nB(ii) = length(find(im.edges(:,1)==ii | im.edges(:,2)==ii));
end

%% Begin centering

if ~visualize_flag
    hwait = waitbar(0,'Centering XYZ...');
end

for jNode=1:length(nLst) %1:nNodes
    % Waitbar
    if ~visualize_flag
        waitbar(jNode/length(nLst),hwait);
    end
    
    % 
    iNode = nLst(jNode);
    
    % If the segment contains b/w zero or one bifurcations
    if nB(iNode)<=2
        [lstE,lstC] = find(im.edges==iNode);

        pos1 = im.nodes(iNode,:);
        pos2 = im.nodes(im.edges(lstE(1),mod(lstC(1),2)+1),:);
        if length(lstE)==2
            pos0 = im.nodes(im.edges(lstE(2),mod(lstC(2),2)+1),:);
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
                    Isub(iy,ix) = angio( iiy, iix, iiz );
                end
            end

            Isub = Isub .* (Isub >= im_thresh); % angio could add a continuity condition here
            % using imfill
            IsubSum = max(sum(Isub(:)),1)+eps;
            [xx,yy] = meshgrid(xLst,yLst);
            posM(1) = sum(xx(:).*Isub(:))/IsubSum;
            posM(2) = sum(yy(:).*Isub(:))/IsubSum;

            if visualize_flag
                figure(1)
                subplot(2,1,1)
                imagesc(angio(:,:,round(pos1(3))))
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
            nodes(iNode,:) = [iix iiy iiz];
            %        pause(0.1)
        end % r>0
    end  % end if nB(iNode)<=2
end % end loop on iNode

im.nodes = nodes;

if ~visualize_flag
    close(hwait)
end
