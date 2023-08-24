function [pMin,pMax,Tree] = plotTreeSeg( E, V, pt, num, flag, im, pMin, pMax, Tree, edgeSegN )

if num==0
    return
end

if ~exist('edgeSegN')
    edgeSegN = ones(size(E,1),1);
end

if flag(1)
    hWait = waitbar(0,'Generating Tree...');
    
    figure(10);
    clf
    
    hold off
    plot3(V(pt,2),V(pt,1),V(pt,3),'r*')
    hold on
    
    pMin = [1000 1000 1000];
    pMax = [1 1 1];
    Tree = [];
end

lst = find( E(:,1)==pt | E(:,2)==pt );

for ii=1:length(lst)
    jj = lst(ii);
    hl=plot3( V(E(jj,:),2), V(E(jj,:),1), V(E(jj,:),3), '.-' );
    set(hl,'ButtonDownFcn',sprintf('plotTreeSeg_select(%d)',jj ));
    if edgeSegN(jj)==0
        set(hl,'color','g');
    elseif edgeSegN(jj)==2
        set(hl,'color','r');
    end
    Tree(end+1) = jj;
    
    for kk=1:3
        for ll=1:2
            if V(E(jj,ll),kk)>pMax(kk)
                pMax(kk) = V(E(jj,ll),kk);
            end
            if V(E(jj,ll),kk)<pMin(kk)
                pMin(kk) = V(E(jj,ll),kk);
            end
        end
    end
    
    En = E;
    En(jj,:) = 0;

    v2 = setdiff(E(jj,:),pt);    
    lstE = find(En(:,1)==v2 | En(:,2)==v2);
    nConn1 = setdiff(unique(En(lstE,:)),v2);
    
    [pMin,pMax,Tree] = plotTreeSeg( En, V, v2, num-1, 0, [], pMin, pMax, Tree, edgeSegN );
    
end

if flag(1)
    close(hWait)
    
    if length(flag)<2
        thresh = 1;
    else
        thresh = flag(2);
    end

    pMin = max(floor(pMin)-20,1);
    pMax = ceil(pMax)+20;
    pMax(1) = min(pMax(1),size(im,2));
    pMax(2) = min(pMax(2),size(im,1));
    pMax(3) = min(pMax(3),size(im,3));
    if exist('im')
        [f,v] = isosurface(im(pMin(2):pMax(2),pMin(1):pMax(1),pMin(3):pMax(3)),thresh   );
        v2 = v+ones(size(v,1),1)*(pMin-1);
        v2 = v2(:,[2 1 3]);
        h = patch('Faces',f,'Vertices',v2);
        set(h,'linestyle','none')
        set(h,'facealpha',.15)
        set(h,'hittest','off');
        axis image
        rotate3d on

%         figure(11)
%         imagesc( max(im(pMin(2):pMax(2),pMin(1):pMax(1),pMin(3):pMax(3)),[],3), [0 1])
%         axis image
%         colormap(gray)
%         hold on
%         plot(V(pt,1)-pMin(1)+1,V(pt,2)-pMin(2)+1,'r*')
% 
%         lst = find( E(:,1)==pt | E(:,2)==pt );
%         for ii=1:length(lst)
%             jj = lst(ii);
%             plot( round(V(E(jj,:),1)-pMin(1)+1), round(V(E(jj,:),2)-pMin(2)+1), 'r.-' )
%         end
%         hold off
% 
%         for jj=1:length(lst)
%             x = repmat([-10:10],[21 1]);
%             y = repmat([-10:10]',[1 21]);
%             z = zeros(21,21);
%             xx = x(:);
%             yy = y(:);
%             zz = z(:);
% 
%             r1 = V(E(lst(jj),1),:);
%             r2 = V(E(lst(jj),2),:);
%             dr = r1-r2;
%             dx = r1(1) - r2(1);
%             dz = r1(3) - r2(3);
%             r = norm(dr);
%             rho = norm(dr(1:2));
%             beta = acos(dz/r);
%             gamma = -acos(dx/rho);
%             Ry = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
%             Rz = [cos(gamma) sin(gamma) 0; -sin(gamma) cos(gamma) 0; 0 0 1];
%             foo = Ry*[xx';yy';zz'];
%             foo = Rz*foo;
%             xx = foo(1,:)';
%             yy = foo(2,:)';
%             zz = foo(3,:)';
%             [ny nx nz] = size(im);
%             %    lst = round((yy+V(pt,2)) + (xx+V(pt,1))*ny + (zz+V(pt,3))*nx*ny);
%             for ii=1:441
%                 A(ii) = im(round(yy(ii)+V(pt,2)),round((xx(ii)+V(pt,1))),round((zz(ii)+V(pt,3))));
%             end
%             figure(12);
%             subplot(1,length(lst),jj)
%             imagesc(reshape(A,[21 21]))
%             axis image
%             hold on
%             plot(11,11,'r*')
%             hold off
%             colormap(gray)
%             title( sprintf('b=%.2f   g=%.2f',beta*180/pi,gamma*180/pi) )
%        end
    end
end    