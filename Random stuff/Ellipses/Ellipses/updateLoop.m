lst = find(im2.nB==1);
lstS = im2.nodeSegN(lst);
%lstDS = find(im2.segLen(lst)<10); % list dangling segments
lstDS = find(im2.segLen(lstS)<im2.segDiam(lstS));
lstDS = lstS(lstDS);

lstDN=(find(ismember(im2.nodeSegN,lstDS)==1)); % list dangling nodes
newPos = zeros(length(lstDN),3);
newPosDel = zeros(length(lstDN),3);
for ii=1:length(lstDN)
    [s,valid] = EllipseFit3DConstrained_dab( I/32,im2.nodePos(lstDN(ii),1),im2.nodePos(lstDN(ii),2),im2.nodePos(lstDN(ii),3) );
    foo=[s.mu-im2.nodePos(lstDN(ii),:)']';
    disp(sprintf('Node %d of %d: updated %.1f   (%.1f %.1f %.1f)', ii, length(lstDN), norm(foo), foo(1), foo(2), foo(3) ))
    newPos(ii,:) = s.mu';
    newPosDel(ii,:) = foo;
end


% Update im2.nodePos
% but don't let update be greater than segment radius
for ii=1:length(lstDN)
    rad = im2.segDiam(im2.nodeSegN(lstDN(ii))) / 2;
    dstep = norm(newPosDel(ii,:));
    if dstep>rad
        dstep2 = rad;
    else
        dstep2 = dstep;
    end
    im2.nodePos(lstDN(ii),:) = im2.nodePos(lstDN(ii),:) + newPosDel(ii,:)*dstep2/dstep;
end


im2.gui.pushbuttonRegraphNodes = 'on';

