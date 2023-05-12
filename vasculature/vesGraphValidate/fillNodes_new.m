function [nodePos, nodeEdges,validatedNodes,verifiedEdges] =...
    fillNodes_new(nodePos, nodeEdges,validatedNodes,verifiedEdges,vox_xy, vox_z)
%%% After downsampling the nodes, 

% Voxel dimensions
hxy = vox_xy;
hz = vox_z;

nN = size(nodePos,1);
nE = size(nodeEdges,1);
hwait = waitbar(0,'Filling in sparse edges');
nEo = nE;
for iE = 1:nEo
    if ~rem(iE,1000)
        waitbar(iE/nEo,hwait)
    end
    
    n1 = nodeEdges(iE,1);
    n2 = nodeEdges(iE,2);
    
    if validatedNodes(n1)==0 && validatedNodes(n2)==0
        len = abs( nodePos(n1,:) - nodePos(n2,:) );
        if sum(len>3*[hxy hxy hz])
            nStep = max(floor(len ./ [hxy hxy hz]) - 1);
            rStep = (nodePos(n2,:) - nodePos(n1,:)) / (nStep+1);
            pos = nodePos(n1,:);
            for jj=1:nStep
                pos = pos + rStep;
                nodePos(end+1,:) = pos;
                validatedNodes(end+1) = 0;
                nN = nN + 1;
                nodeEdges(end+1,:) = [n1 nN];
                verifiedEdges(end+1) = 0;
                nE = nE + 1;
                n1 = nN;
            end
            nodeEdges(iE,:) = [nN n2];
        end
    end
end
close(hwait)