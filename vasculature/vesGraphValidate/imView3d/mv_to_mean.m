function im = mv_to_mean(im, Ithresh)
% Move towards mean of neighboring nodes
nodes = im.nodes;
nLst = 1:size(nodes,1);
hwait = waitbar(0,'Moving nodes towards mean of neighboring nodes');
for iii= 1:length(nLst) %1:nN
    ii = nLst(iii);
    waitbar(iii/length(nLst),hwait)
    eLst = find(im.edges(:,1)==ii | im.edges(:,2)==ii);
    nLst2 = setdiff(unique(im.edges(eLst,:)), ii);
    if length(nLst2)>1
        pos0 = max(im.nodes(ii,:),1);
        posC = mean(im.nodes(nLst2,:),1);
        posN = pos0 + (posC-pos0) / max(norm(posC-pos0),1);
%            if im.angio(round(pos0(2)),round(pos0(1)),round(pos0(3)))>=Ithresh
            if im.angio(round(posN(2)),round(posN(1)),round(posN(3)))>=Ithresh
                nodes(ii,:) = posN;
            end
%            else
%                nodes(ii,:) = posN;
%            end
    end
end
im.nodes = nodes;
close(hwait)
end