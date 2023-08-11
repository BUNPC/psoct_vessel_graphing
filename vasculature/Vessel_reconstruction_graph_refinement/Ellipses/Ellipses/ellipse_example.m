%% Example code of marching ellipse
%{
This is an outline from Sava.
%}

load matlab;

%% 
s = EllipseFit3DConstrained_dab( V,72,77,40,1);

%% 

ShowLocalDataWithSE(V,s,'foo');

%% 

% how to find points close to ellipsoid center, inside a cylinder defined
% by ellipsoid radius and height equal to fill factor length

% fill factor = f (it should be hard coded somewhere in imView3d - length
% between nodes

% s.a2 is radius in pixels of ellipsoid longest dimension (along vessel)
% s.a1 and s.a3 are other two ellipsoid radiuses

% Following code will generate (XYZ)radVec matrices. Each contains two point
% coordinates - one from ellipsoid center (s.mu) and other at the end of
% ellipse's X, Y, or Z radius

XradiusVec = [[0 s.a1]' [0 0]' [0 0]']';
YradiusVec = [[0 0]' [0 s.a2]' [0 0]']';
ZradiusVec = [[0 0]' [0 0]' [0 s.a3]']';

R = rotation_quat(s.q);
XradVec = R*XradiusVec + repmat(s.mu,1,size(XradiusVec,2));
YradVec = R*YradiusVec + repmat(s.mu,1,size(YradiusVec,2));
ZradVec = R*ZradiusVec + repmat(s.mu,1,size(ZradiusVec,2));

%% 
% these are just vectors of ellipse radiuses translated to origin:
Xvec = XradVec(:,2)-XradVec(:,1);
Yvec = YradVec(:,2)-YradVec(:,1);
Zvec = ZradVec(:,2)-ZradVec(:,1);

% They can be generated faster than in above code, but lets not think
about this for now...
%%
% for any node at point (X,Y,Z) in vessel stack, we can estimate if it is
% inside cylinder around ellipse's center...

% loop through all nodes in graph...
node = [2, 4, 6] - s.mu; % example signle node (x, y, z) coordinates
are [2, 4, 6]. Needed to subtract s.mu to get node in respect to
ellipse center

if abs(node*Yvec)<f && ( abs( node*Xvec + node*Zvec )< max(s.a1, s.a3)
), % looks at scalar vector products to decide if point is inside
cylinder...
   % point is inside cylinder, move the node and set its flag
   nodeflag = 1;
else
   % point is outside cylinder
   % don't do anything
end;


You will need the last part of it - after calculation of 's'. Please
let me know how it is going... One thing that is wrong here is that
ellipse seems to be generated without respect to physical dymensions
of voxels. I didn't try to correct this, so numbers are all in pixels.
Lets first see how it works like this, and then we can see if we need
to improve it further...

Best regards
Sava
