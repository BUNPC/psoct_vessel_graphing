function [II, info, trkPos] = DTItrack( I, p1, cxyz, nTracksMax, biDirectional )

% can I make this better at turning corners?

% make sure max(I)<=1
% I2 = im.Istk;
% I2 = min(I2/2,1);

% need to filter BG
% need gaussian filter
% I = BGFilter( I2, 1.2);
% I = gaussianFilter( I, 1 );

% need a viewer of tracks overlayed on image
% imView3d(I,II)

%    x=154; y=15; z=16; x1o=0;y1o=1;z1o=0; % tile 33
%    x=212; y=170; z=24; x1o=-1;y1o=0;z1o=0;
% xi=30; yi=60; zi=18; xio=1;yio=1;zio=0;
%    x=120;y=240;z=20; x1o=0;y1o=-1;z1o=0;
% xi=140;yi=180;zi=50;xio=1;yio=-1;zio=0;  % Seed 1
% xi=207;yi=128;zi=39; xio=0;yio=0;zio=0; % seed 2

xi = p1(1);
yi = p1(2);
zi = p1(3);

xio = cxyz(1);
yio = cxyz(2);
zio = cxyz(3);



Imin = 0.1;
Imax = 1;

f(1) = Imin;
f(2) = 0.5;

h = 5;
hz = 2;
stepSize = 1;
cPersistence = 0.8;
nStepMax = 10000;
%nTracksMax = 100;

info.Imin = Imin;
info.Imax = Imax;
info.f = f;
info.h = h;
info.hz = hz;
info.stepSize = stepSize;
info.cPersistence = cPersistence;
info.nStepMax = nStepMax;
info.nTracksMax = nTracksMax;
info.posO = [xi yi zi];
info.cO = [xio yio zio];

hh = 2*h+1;
hhz = 2*hz+1;
II = zeros(size(I));

%figure(1)
%imagesc(I,[Imin Imax])
%colormap(gray)

nx1 = h+1;
ny1 = h+1;
nz1 = hz+1;

nx2 = size(I,2)-h-1;
ny2 = size(I,1)-h-1;
nz2 = size(I,3)-hz-1;

trkPos = zeros( nStepMax, 3, nTracksMax );

for nTracks = 1:nTracksMax
    x = xi; y = yi; z = zi;
    sgn = sign(biDirectional^nTracks+eps);
    x1o = xio*sgn; y1o = yio*sgn; z1o = zio*sgn;
    
    xx = round(x);
    yy = round(y);
    zz = round(z);

    IItemp = zeros(size(I));
    
    nStep = 1;
    nStep2 = 1;
    rnd = rand(nStepMax,1);
    while nStep<nStepMax & xx>nx1 & xx<nx2 & yy>ny1 & yy<ny2 & zz>nz1 & zz<nz2
        nStep = nStep + 1;

        p = (max(I(yy-h:yy+h,xx-h:xx+h,zz-hz:zz+hz) - f(1),0)).^f(2);

        sp = sum(p(:));
        if sp>0
            IItemp(yy,xx,zz) = IItemp(yy,xx,zz) | 1;

            [p,lst] = inNormSortP( p, sp );
%            p = p(:)/sp;
%            [p,lst] = sort(p);

            ii = floor(length(p)/2);
            h2 = floor(ii/2);
            rndNum = rnd(nStep);
            ii = inFindIdx( p, ii, h2, rndNum );
%            while h2>1
%                if sum(p(1:ii))>rnd(nStep)
%                    ii = ii - h2;
%                else
%                    ii = ii + h2;
%                end
%                h2 = ceil(h2 / 2);
%            end
%            ii = min(max(ii,1),length(lst));

            x1 = ceil(mod(lst(ii),hh*hh)/hh) - h - 1;
            y1 = mod(lst(ii),hh) - h - 1;
            z1 = ceil(lst(ii)/(hh*hh)) - hz - 1;
            r = norm([y1 x1 z1])+eps;
            x1 = x1/r;
            y1 = y1/r;
            z1 = z1/r;
            dcor = sign( [x1 y1 z1]*[x1o y1o z1o]' );
            if dcor>=0
                x = x + stepSize*x1;
                y = y + stepSize*y1;
                z = z + stepSize*z1;
                xx = round(x);
                yy = round(y);
                zz = round(z);

                trkPos( nStep2, :, nTracks ) = [xx yy zz];

                x1o = x1o + cPersistence*x1;
                y1o = y1o + cPersistence*y1;
                z1o = z1o + cPersistence*z1;
                r = norm([x1o y1o z1o]);
                x1o = x1o / r;
                y1o = y1o / r;
                z1o = z1o / r;
                
                nStep2 = nStep2 + 1;
            end
        else
            x = x - stepSize*x1;
            y = y - stepSize*y1;
            z = z - stepSize*z1;
            xx = round(x);
            yy = round(y);
            zz = round(z);
            
            nStep = nStep - 1;

%            II(yy,xx,zz) = II(yy,xx,zz) - 1;
        end
        

    end
    
    II = II + IItemp;
    
%     if mod(nTracks,10)==0
%         III = zeros(size(I));
%         lst = find(I>Imin);
%         III(lst) = 1;
%         III = III + II;
% 
%         figure(2)
%         imagesc(III,[0 4])
%     end
end


%Imin2 = Imin;
%Imin2 = Imax * 0.4;
%III = zeros(size(I));
%lst = find(I>Imin2);
%III(lst) = 1;
%III = III + II;


return






function [p,lst] = inNormSortP( p, sp )
p = p(:)/sp;
[p,lst] = sort(p);



function ii = inFindIdx( p, ii, h2, rndNum )
while h2>1
    if sum(p(1:ii))>rndNum
        ii = ii - h2;
    else
        ii = ii + h2;
    end
    h2 = ceil(h2 / 2);
end
ii = min(max(ii,1),length(p)); % length(lst)
