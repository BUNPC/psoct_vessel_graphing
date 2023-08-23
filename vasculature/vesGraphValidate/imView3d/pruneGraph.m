function [nodeEdges,nB2] = pruneGraph( im );

% repeat this pruning at the 2nd level until there is no change
% this leaves orphan nodes that need to be pruned.

nNodes = size(im.nodePos,1);

nE0 = size(im.nodeEdges,1);

nB=zeros(1,nNodes);
for ii=1:nNodes
    nB(ii)=length(find(im.nodeEdges(:,1)==ii | im.nodeEdges(:,2)==ii)); 
end
nBtmp = nB;

lst = find(nB>2);

edgeFlag = ones(1,size(im.nodeEdges,1));

hWait = waitbar(0,'Pruning edges...');
for iiNidx = 1:length(lst)
    nIdx = lst(iiNidx);
    waitbar(find(nIdx==lst)/length(lst),hWait);
    lstE = find(im.nodeEdges(:,1)==nIdx | im.nodeEdges(:,2)==nIdx);
    nConn1 = setdiff(unique(im.nodeEdges(lstE,:)),nIdx);
    flagBreak = 0;
    for ii=1:length(nConn1)
        lstE1 = find(im.nodeEdges(:,1)==nConn1(ii) | im.nodeEdges(:,2)==nConn1(ii));
        nConn2 = setdiff(unique(im.nodeEdges(lstE1,:)),[nIdx nConn1(ii)]);


        for jj=1:length(nConn2)
            lstE2 = find(im.nodeEdges(:,1)==nConn2(jj) | im.nodeEdges(:,2)==nConn2(jj));
            nConn3 = setdiff(unique(im.nodeEdges(lstE2,:)),[nConn1(ii) nConn2(jj)]);

            %        [nIdx nConn1(ii) nConn2(jj) nConn3]
            if ~isempty(find(nConn3==nIdx))
                lst2 = find(im.nodeEdges(lstE2,1)==nIdx | im.nodeEdges(lstE2,2)==nIdx);
                if length(lst2)>1
                    warning( 'DAB LOOK AT THIS!!!!' )
                end
                kk = 1;
                if nBtmp(nConn2(jj))>2
                    edgeFlag(lstE2(lst2(kk))) = 0;
                    nBtmp(im.nodeEdges(lstE2(lst2(kk)),1)) = nBtmp(im.nodeEdges(lstE2(lst2(kk)),1)) -1;
                    nBtmp(im.nodeEdges(lstE2(lst2(kk)),2)) = nBtmp(im.nodeEdges(lstE2(lst2(kk)),2)) -1;

                    flagBreak = 1;
%                    iiNidx = iiNidx - 1;
                elseif nBtmp(nConn1(ii))==2
                    % I had added this for a reason, but it was bad.
                    % I wish I could remember why I added it
                    % I remember, It was to deal with 3 edge loops at end
                    % points.. i.e. consider edges = [0 1;1 2; 2 3; 3 4;4 2]
                    % This is a kite with a loop at the end. We don't want
                    % that loop at the end
%                     edgeFlag(lstE2(lst2(kk))) = 0;
%                     edgeFlag(lstE1) = 0;
%                     nBtmp(nIdx) = nBtmp(nIdx) - 2;
% 
%                     flagBreak = 1;
%                    iiNidx = iiNidx - 1;
                end
%                 for kk=1:length(lst2)
%                     if nBtmp(im.nodeEdges(lstE2(lst2(kk)),1))>2 & nBtmp(im.nodeEdges(lstE2(lst2(kk)),2))>2
%                         edgeFlag(lstE2(lst2(kk))) = 0;
%                         nBtmp(im.nodeEdges(lstE2(lst2(kk)),1)) = nBtmp(im.nodeEdges(lstE2(lst2(kk)),1)) -1;
%                         nBtmp(im.nodeEdges(lstE2(lst2(kk)),2)) = nBtmp(im.nodeEdges(lstE2(lst2(kk)),2)) -1;
% 
%                         flagBreak = 1;
%                         iiNidx = iiNidx - 1;
%                     end
%                     if flagBreak
%                         break
%                     end
%                 end
            end
            if flagBreak
                break
            end
        end
        
        if flagBreak
            break
        end
    end
end

nodeEdges = im.nodeEdges(find(edgeFlag==1),:);

%nodeFlag = zeros(1,size(nodeEdges,1));


nB2=zeros(1,nNodes);
for ii=1:nNodes
    nB2(ii)=length(find(nodeEdges(:,1)==ii | nodeEdges(:,2)==ii)); 
end

close(hWait)

figure
hist([nB' nB2'],[0:20])
title( sprintf('%d edges reduced to %d',nE0,size(nodeEdges,1)) )
    
