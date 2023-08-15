function [s, valid] = EllipseFit3DConstrained(Vol, x, y, z, ax, ay, az)
% Contrained ellipsefit, where the solution is found on the plane plane
% perpendicular to the vector [ax, ay, az].
% Vol is the volume data, and [x, y, z] are initial position vector.
% the output is a structure s containing various ellipse parameters.
% Note that I have constrained the ellipsoid to 8 degrees of freedom that
% are 3 position , 3 orientation and 2 scales. The major axis (along the
% vessel is given by AS_RATIO*max( minor axes). 
% Author: Amit Mukherjee

global PARAM;
[dim.y, dim.x, dim.z] = size(Vol);

s.a1 = 4;
s.a2 = 4;
s.a3 = 6;
s.q = [0 1 0 0]';
s.mu = [x y z]';
s.L = 0;
s.MAD = 0;
ConstAxis = [ax, ay, az]';

Iteration  = 40;
valid = 1;

%initialise parameters here
PARAM.dt_u = 500*[1,1,1];
PARAM.dt_a = 100*[1,1,1]';
PARAM.dt_q = 500*[1 1 1 1]';
PARAM.sign_Q = [0,0,0,0]';
PARAM.sign_A = [0,0,0]';
PARAM.sign_U = [0,0,0];

t = SphericalTesselation(29);
flag = 0;
i =0;
j=0;
%for i = 1:Iteration
while flag<2
    i = i +1;
    j=j+1;
    sLo = s.L;
     %if mod(i,5) == 0
        s = UpdateIntensityMedian(Vol,s);
     %end
    [t,s] = getSurface(s,t);
    t = UpdateForcesChanVese(Vol,s,t);
    s = UpdateOrientation(s,t);
    s = UpdateScale(s,t);
    if flag>0
        R = rotation_quat(s.q);
        ConstAxis = R(:,3);
        s = UpdatePosition(s,t, ConstAxis);
    end
    if abs(s.L-sLo)/sLo < 1e-7 && j>10
        flag = flag + 1;
        j=0;
    end
    fprintf('Iter: %d (%.3f)  Rate = [%.02f, %.02f, %.02f, %.02f, %.02f, %.02f %.02f, %.02f, %.02f]\n', i, s.L,...
        PARAM.dt_u(1), PARAM.dt_u(2), PARAM.dt_u(3), PARAM.dt_a(1), PARAM.dt_a(2), ...
        PARAM.dt_q(1), PARAM.dt_q(2), PARAM.dt_q(3), PARAM.dt_q(4) );
    ShowLocalDataWithSE(Vol,s, ['Iter: ',num2str(i), ' Likelihood: ',num2str(s.L)]); 
    drawnow
%    pause(0.5);

%   if ExitCriteria(s,t, PARAM, i)
%         break;
%    end
end
s = UpdateIntensityMedian(Vol,s);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = UpdatePosition(s,t,ca)
global PARAM

delta_u(1) = -PARAM.dt_u(1)*(t.normal(1,:)*t.Lhd')/t.numFaces;
delta_u(2) = -PARAM.dt_u(2)*(t.normal(2,:)*t.Lhd')/t.numFaces;
delta_u(3) = -PARAM.dt_u(3)*(t.normal(3,:)*t.Lhd')/t.numFaces;

% remove the component of displacement along constraint axis
delta_u = delta_u - ca'*(delta_u*ca);

% limit the maximum value of update
delta_u = sign(delta_u) .* min(abs(delta_u), 2);
%fprintf('dmu = [%.02f, %.02f, %.02f]\n',delta_u(1), delta_u(2), delta_u(3));

s.mu(1) = s.mu(1) + delta_u(1);
s.mu(2) = s.mu(2) + delta_u(2);
s.mu(3) = s.mu(3) + delta_u(3);

PARAM.dt_u = (1-0.2*double(sign(delta_u)~=PARAM.sign_U)).*PARAM.dt_u;
PARAM.sign_U = sign(delta_u);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = UpdateScale(s,t)
% da = R*(s/a) * n * Lhd * area
global PARAM 
MAX_VESSEL_WIDTH = 20;
MIN_VESSEL_WIDTH = 2;
AS_RATIO = 1.5;

R = rotation_quat(s.q);
da1 = 0;
da2 = 0;
for i = 1:t.numFaces
    da1 = da1 + ( R(:,1)*t.local_centroid(1,i)/s.a1 )' * t.normal(:,i) * t.Lhd(i);
    da2 = da2 + ( R(:,2)*t.local_centroid(2,i)/s.a2 )' * t.normal(:,i) * t.Lhd(i);
end
%%%%%%%%%%%%%%%%%
da1 = da1/t.numFaces * PARAM.dt_a(1);
da2 = da2/t.numFaces * PARAM.dt_a(2);

%limit the max absolute value of da1/da2 between -2 and 2 pixels
if (da1 < -1.0), da1 = -1.0; end
if (da1 >  1.0), da1 =  1.0; end
if (da2 < -1.0), da2 = -1.0; end
if (da2 >  1.0), da2 =  1.0; end
fprintf(' a = [%.02f, %.02f]', s.a1,  s.a2);

%multiply by dt_a
s.a1 = s.a1 - da1 ;
s.a2 = s.a2 - da2 ;

%limit the min-max value of a1/a2 between MIN_VESSEL_WIDTH  and
%MAX_VESSEL_WIDTH and their ratio < 2

if s.a1 > MAX_VESSEL_WIDTH, s.a1 = MAX_VESSEL_WIDTH; end
if s.a1 < MIN_VESSEL_WIDTH, s.a1 = MIN_VESSEL_WIDTH; end
if s.a2 > MAX_VESSEL_WIDTH, s.a2 = MAX_VESSEL_WIDTH; end
if s.a2 < MIN_VESSEL_WIDTH, s.a2 = MIN_VESSEL_WIDTH; end

if (s.a1 / (s.a2+0.0001)) > 2.0,    s.a1 = s.a2 * 2.0; end
if (s.a2 / (s.a1+0.0001)) > 2.0,    s.a2 = s.a1 * 2.0; end

% approximate a3
s.a3 = AS_RATIO * max(s.a1, s.a2);

if(da1 * PARAM.sign_A(1) <= 0 ), PARAM.dt_a(1) = PARAM.dt_a(1)*0.8; end
if(da2 * PARAM.sign_A(2) <= 0 ), PARAM.dt_a(2) = PARAM.dt_a(2)*0.8; end

PARAM.sign_A(1) = da1;
PARAM.sign_A(2) = da2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = UpdateOrientation(s,t)

global PARAM

R = rotation_quat(s.q);

s0 = s.q(1);
v1 = s.q(2);
v2 = s.q(3);
v3 = s.q(4);
G = 2*[[-v1 s0 v3 -v2];[-v2 -v3 s0 v1];[-v3 v2 -v1 s0]];
 
delta_q = zeros(4,1);
for i = 1:t.numFaces
    p1 = t.local_centroid(1,i); 
    p2 = t.local_centroid(2,i); 
    p3 = t.local_centroid(3,i); 
    p = [0 -p3  p2; ...
        p3  0  -p1; ...
        -p2 p1  0];
    C = R*p*G;
    delta_q = delta_q + C'*t.normal(:,i)*t.Lhd(i);
end

delta_q = delta_q ./t.numFaces;
delta_q =  delta_q .* PARAM.dt_q;

delta_q = sign(delta_q).*min(abs(delta_q),0.05);
s.q = s.q + delta_q;
s.q = s.q/norm(s.q);

if sum(sign(delta_q)~=PARAM.sign_Q)
    PARAM.dt_q = PARAM.dt_q*.8;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = ExitCriteria(s, t, PARAM, iter)
res = 0;
if max(s.a1,s.a2) > PARAM.max_width
    res = 1;
end

if PARAM.toggle_count_A > 20 || sum(PARAM.net_change_A) < 0.1
    res = 1;
end

if PARAM.toggle_count_Q > 50 && abs(sum(PARAM.net_change_Q(:))) < 0.01
    res = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = UpdateIntensityMedian(Vol,s)
[bg_data, fg_data, n_bg, n_fg] = GetInsideOutsideVoxels(Vol, s);
F_est = median(fg_data);
B_est = median(bg_data);
% F_est = mean(fg_data);
% B_est = mean(bg_data);

L = F_est - B_est;

if L > s.L
    s.f = F_est;
    s.b = B_est;
    s.L = L;
end
if L < 0.001
    mn = 0.5 * (F_est + B_est);
    s.f = mn + 0.001;
    s.b = mn - 0.001;
    s.L = 0.001;
end
s.b = max( s.b , 0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bg_data, fg_data, bg_ndx, fg_ndx] = GetInsideOutsideVoxels(Vol, Segment)
[d1,d2,d3] = size(Vol);
A = Segment.a3;

lx = round(Segment.mu(1) - A); hx = round(Segment.mu(1) + A);
ly = round(Segment.mu(2) - A); hy = round(Segment.mu(2) + A);
lz = round(Segment.mu(3) - A); hz = round(Segment.mu(3) + A);

R = rotation_quat(Segment.q);
% R = R'; %inverse

bg_data = zeros(1,ceil(9*A*A*A/2)); bg_ndx = 0;
fg_data = zeros(1,ceil(9*A*A*A/2)); fg_ndx = 0;
InBoundary = 0;
OutBoundary = 0;

for x = lx : hx
    for y = ly : hy
        for z = lz : hz
            X = [x, y, z]';
            xx = (R(1,:)* (X - Segment.mu))/Segment.a1;
            yy = (R(2,:)* (X - Segment.mu))/Segment.a2;
            zz = (R(3,:)* (X - Segment.mu))/Segment.a3;
            if ( (x>0)&&(x<=d2) && (y>0)&&(y<=d1) && (z>0)&&(z<=d3) )
                InBoundary = InBoundary + 1;
                %v(  y-ly+1 ,x-lx+1 , z-lz+1) = Vol(y,x,z);
                D = xx.^2 +yy^2 + zz^2;
                if D > 1
                    bg_ndx = bg_ndx + 1;
                    bg_data(bg_ndx) = Vol(y,x,z);
                else
                    fg_ndx = fg_ndx + 1;
                    fg_data(fg_ndx) = Vol(y,x,z);
                end
            else
                    OutBoundary = OutBoundary + 1;
            end
        end
    end
end

fg_data(fg_ndx+1:end) = [];
bg_data(bg_ndx+1:end) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,s] = getSurface(s,t)
t.scaled = diag([s.a1,s.a2,s.a3])*t.vertex;
%vertex scaled, rotrated and translated in global coordinates
R = rotation_quat(s.q);
t.global = R*t.scaled + repmat(s.mu,1,size(t.vertex,2));

for k = 1:t.numFaces
    %local centroid from unrotated mesh
    tmp1 = [t.scaled(1,t.n1(k));t.scaled(2,t.n1(k));t.scaled(3,t.n1(k))]; %node 1
    tmp2 = [t.scaled(1,t.n2(k));t.scaled(2,t.n2(k));t.scaled(3,t.n2(k))]; %node 2
    tmp3 = [t.scaled(1,t.n3(k));t.scaled(2,t.n3(k));t.scaled(3,t.n3(k))]; %node 3
    t.local_centroid(:,k) = (tmp1+tmp2+tmp3)/3;

    %global centroid     
    tmp1 = [t.global(1,t.n1(k)), t.global(2,t.n1(k)), t.global(3,t.n1(k))]; %global 1
    tmp2 = [t.global(1,t.n2(k)), t.global(2,t.n2(k)), t.global(3,t.n2(k))]; %global 2
    tmp3 = [t.global(1,t.n3(k)), t.global(2,t.n3(k)), t.global(3,t.n3(k))]; %global 3

    t.global_centroid(:,k) = (tmp1+tmp2+tmp3)/3;
    dif1 = tmp1 - tmp2; %sides of triangle DONT change this
    dif2 = tmp3 - tmp1;
    %t.normal = cross(dif1',dif2',1);  %must be 3 in rows n in columns
    t.normal(1,k) = ( dif1(2)*dif2(3) - dif1(3)*dif2(2) );
    t.normal(2,k) = (-dif1(1)*dif2(3) + dif1(3)*dif2(1) );
    t.normal(3,k) = ( dif1(1)*dif2(2) - dif1(2)*dif2(1) );

    mag = sqrt(sum(t.normal(:,k).^2));
    t.area(k) = 0.5*mag;
    t.normal(:,k) = t.normal(:,k)/mag;
end
