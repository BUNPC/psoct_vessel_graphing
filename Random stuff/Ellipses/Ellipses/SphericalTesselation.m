function t = SphericalTesselation(N)
%this is tested working 
TEST = 0;
nu = (N-1)/2;
nv = N+1;
t.numVertex = nu*nv + 2;
t.numFaces = nv*(nu-1)*2 + 2*nv;


step = 2*pi/(N+1);
%for the new tesselation method
%for u leave the end points. and v make the complete circle
%[u,v]=meshgrid( -pi/2+step: step : pi/2-step, -pi : step : pi-step );
q = 0;
u = zeros(1,nu*nv);
v = zeros(1,nu*nv);
for iu = 0:nu-1
    for iv = 0:nv-1
        q = q+1;
        u(q) = step*(iu+1) - pi/2;
        v(q) = step*iv - pi;
    end
end
% [u(:)-u1(:) v(:)-v1(:)]
%add the end points
t.u = [-pi/2; u(:); pi/2 ]';    
t.v = [-pi ; v(:); pi ]';       

t.Cu = cos(t.u);
t.Cv = cos(t.v);
t.Su = sin(t.u);
t.Sv = sin(t.v);
t.vertex = [t.Cu'.*t.Cv', t.Cu'.*t.Sv' , t.Su'];


%faces (tringles computation)
t.n1 = zeros(1,t.numFaces);
t.n2 = zeros(1,t.numFaces);
t.n3 = zeros(1,t.numFaces);

ndx = 0;
%both the triangles are counter clockwise
%IMPORTANT: NORMALS SHOULD POINT OUTWARDS
for i = 1:nv*(nu-1)
    if mod(i,nv) 
        %regular
        A = i + 1       ;
        B = i+nv + 1   ;
        C = i+nv+1 + 1 ;
        D = i+1 + 1    ;
    else
        %wrapping
        A = i + 1      ;
        B = i+nv + 1   ;
        D = i-nv+1 + 1 ;
        C = i+1 + 1    ;
    end
%   first triangle ABD
    ndx = ndx+1;
    t.n1(ndx) = A;
    t.n2(ndx) = B;
    t.n3(ndx) = D;
    
%   second triangle DBC
    ndx = ndx+1;
    t.n1(ndx) = D;
    t.n2(ndx) = B;
    t.n3(ndx) = C;
    
    
%    disp(([t.n1(ndx-1) t.n2(ndx-1) t.n3(ndx-1);  t.n1(ndx) t.n2(ndx) t.n3(ndx)]));
%     disp(t.v([t.n1(2*i-1) t.n2(2*i-1) t.n3(2*i-1);  t.n1(2*i) t.n2(2*i) t.n3(2*i)])*180/pi);
%     disp(t.u([t.n1(2*i-1) t.n2(2*i-1) t.n3(2*i-1);  t.n1(2*i) t.n2(2*i) t.n3(2*i)])*180/pi);
end

% %offset of one due to leaving the first point
% t.n1 = t.n1+1;
% t.n2 = t.n2+1;
% t.n3 = t.n3+1;

%test code
if TEST == 1
     f = [t.n1(1:ndx)' t.n2(1:ndx)' t.n3(1:ndx)'] ;
     patch('Faces',f(1:2:end,:),'Vertices',t.vertex, 'facecolor','g');
     patch('Faces',f(2:2:end,:),'Vertices',t.vertex, 'facecolor','r');
end
%end test code
%top pole
for i = [1:nv] + 1
    ndx = ndx+1;
    t.n1(ndx) = 1;
    t.n2(ndx) = i;
    if i ~= nv+1
        t.n3(ndx) = i+1;    
    else
        %wrapping
        t.n3(ndx) = 2;
    end
end
if TEST == 1
    f = [t.n1(ndx-nv+1:ndx)' t.n2(ndx-nv+1:ndx)' t.n3(ndx-nv+1:ndx)'];
    patch('Faces',f(1:end,:),'Vertices',t.vertex, 'facecolor','b');
end
%bottom pole
for i = [1:nv] + nv*(nu-1)+1
    ndx = ndx+1;
    t.n1(ndx) = i;
    t.n2(ndx) = nv*nu + 2;
    if i ~= (nv*nu+1)
        t.n3(ndx) = i+1;    
    else
        %wrapping
        t.n3(ndx) = nv*(nu-1)+2;
    end
end
if TEST == 1
    f = [t.n1(ndx-nv+1:ndx)' t.n2(ndx-nv+1:ndx)' t.n3(ndx-nv+1:ndx)'];
    patch('Faces',f(1:end,:),'Vertices',t.vertex, 'facecolor','c');
    disp(ndx);
    hold on; 
end
%display normals and centroids
if TEST == 1
    for k = 1:t.numFaces
        tmp1 = [t.vertex(t.n1(k),1), t.vertex(t.n1(k),2), t.vertex(t.n1(k),3)]; %vertex 1
        tmp2 = [t.vertex(t.n2(k),1), t.vertex(t.n2(k),2), t.vertex(t.n2(k),3)]; %vertex 2
        tmp3 = [t.vertex(t.n3(k),1), t.vertex(t.n3(k),2), t.vertex(t.n3(k),3)]; %vertex 3

        t.centroid(k,:) = (tmp1+tmp2+tmp3)/3;
        dif1 = tmp1 - tmp2; %sides of triangle DONT change this
        dif2 = tmp3 - tmp1;
        %t.normal = cross(dif1',dif2',1);  %must be 3 in rows n in columns
        t.normal(k,1) = ( dif1(2)*dif2(3) - dif1(3)*dif2(2) );
        t.normal(k,2) = (-dif1(1)*dif2(3) + dif1(3)*dif2(1) );
        t.normal(k,3) = ( dif1(1)*dif2(2) - dif1(2)*dif2(1) );
        
        mag = sqrt(sum(t.normal(k,:).^2));
        t.area(k) = 0.5*mag;
        t.normal(k,:) = t.normal(k,:)/mag;

        plot3(t.centroid(k,1), t.centroid(k,2), t.centroid(k,3), 'r.');
        line(t.centroid(k,1) + [0 t.normal(k,1)], ...
            t.centroid(k,2) + [0 t.normal(k,2)], ...
            t.centroid(k,3) + [0 t.normal(k,3)]);
    end
end
hold off

t.vertex = t.vertex';
