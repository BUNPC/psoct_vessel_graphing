function t = UpdateForcesChanVese(Vol,s,t)
% verify the sign of K later
[dim.y, dim.x, dim.z] = size(Vol);

t.Lhd = ones(1,t.numFaces);
t.F = ones(1,t.numFaces);

for i = 1:t.numFaces
    p = t.global_centroid(:,i);
    l = floor(p);
    h = ceil(p);
    f = p - l;
    if (min(l) <= 0) || (h(1)>dim.x) || (h(2)>dim.y) || (h(3)>dim.z)
        t.Lhd(i) = t.area(i);
        t.F(i) = s.b;
    else
        w000 = Vol(l(2),l(1),l(3));   
        w010 = Vol(h(2),l(1),l(3));   
        w100 = Vol(l(2),h(1),l(3));
        w001 = Vol(l(2),l(1),h(3));   
        w110 = Vol(h(2),h(1),l(3));   
        w101 = Vol(l(2),h(1),h(3));
        w011 = Vol(h(2),l(1),h(3));
        w111 = Vol(h(2),h(1),h(3));
        
        w00 = w000*(1-f(3)) + w001*f(3);
        w01 = w010*(1-f(3)) + w011*f(3);
        w10 = w100*(1-f(3)) + w101*f(3);
        w11 = w110*(1-f(3)) + w111*f(3);
        
        w0 = w00*(1-f(2)) + w01*f(2);
        w1 = w10*(1-f(2)) + w11*f(2);
        
        w  = w0*(1-f(1)) + w1*f(1);
        t.F(i) = w;
        t.Lhd(i) = t.area(i) * (abs(s.f-w) - abs(s.b-w)); 
    end
end
