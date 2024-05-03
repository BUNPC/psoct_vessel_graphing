function I = function_Profile2(x,xdata)
global g

yy = g.yy;
% Ioff = g.Ioff;
% IoffSlope = g.IoffS;


% % Confocal parameters
% load('/autofs/space/vault_003/users/Hui/ScattCoefFitting/TissuePropertiesFitting/ConfocalParameter10x.mat','ConfocalParameters_10x');
% zR = 1.4 * ConfocalParameters_10x(2); % /2.9; % convert back to pixel number
%     % zR = pi w_o^2 / lambda = lambda / (pi NA^2)
%     % wo = lambda / (pi NA)

mus = x(1)/1000; %x(1)/1000;
zf = x(2)+g.zz;
% zf = -6.33 +g.zz;
% zf = -57.48 +g.zz; %48.71, 22.68, -5.33, -6.33, -57.48 
% -57.48, -6.33, -5.33, 22.68, 48.71
%(zoff = 25) 1: -112.15, -36.4; 2: -25.30 (zoff = 15); 3: -65.16; 4: -69.47; 5: -126.09
zR = x(3); % 1: 80.96; 2: 78.96; 3: 72.35; 4: 64.57; 5: 62.63
% zR = 64; %80.96, 78.67, 73, 70, 64
% 62.9614   62.0499   66.6023   73.5360   90.5059

z = xdata; % should I subtract off surface so that z=0 is the surface? DONE

n = 1.3228;
f = 41.35e3; % 4.5e3;  16.67e3; 10x: 41.35e3; 5x: 76.79e3  % um   / 2.9; % 2.9 um/pixel
lambda = 1.3; % um     / 2.9;
%D = 1000; % um    / 2.9; %beam diameter at objective

zp = z - zf; % zero is at the beam waist
pi = 3.14159;

if 1
    % Caroline.
    % I don't know where this comes from
    %     Io = x(4);
    w = ( 1 + ( ( g.zs*3.03 + z - zf ) ./zR ).^2) ;
%     w = ( 1 + ( ( g.zs+z - zf ) ./zR ).^2) ;
    I = ( 1./w ) .* exp(-2*mus*(z+g.zs*3.03));
%     I = ( 1./w ) .* exp(-2*mus*(z));
    %     I = Io * I;
    %     g.Io = Io;
    Io = sum(yy'.*I) / sum(I.^2);
    I = Io.*I;
    g.Io = Io;
    
elseif 0
    % Schmitt
    % This doesn't properly handle the case of air tissue interface as far
    % as I can tell. But should be fine with water immersion and index
    % matching. But as I look more at Schmitt versus Izatt, I realize I
    % can't really derive Schmitt from Izatt. I understand Izatt better and
    % it applies to the general case of having an air tissue interface
    D = x(3);
    I = exp(-2*mus*z) ./ ( 2*(n*zp+f).^2 .* (1 + (n*zp*pi*D^2./(4*lambda*f*(n*zp+f))).^2 ));
    if 0
        Io = sum((yy'-Ioff-IoffSlope*z).*I)/sum(I.^2);
        I = Io.*I + Ioff + IoffSlope*z;
    else
        foo = (IoffSlope/Ioff)*z + 1;
        Io = sum(yy'.*I.*foo-Ioff*I.*foo.^2) / sum(I.^2.*foo.^2);
        I = (Io.*I + Ioff).*foo;
        g.Io = Io;
    end
else
    % Izatt
    % This does properly handle the air tissue interface
    % Says they are derived from Schmitt, but I can't fully see it
    if length(x)==2
        D = 6200; %727;%x(3);
    else
        D = x(3); % Izatt model: dil1x: 6597; dil5x: 6813; dil10x: 6505; 1024, 5xobj, dil5x: 6564
        % 102516: 5xobj: 6321.6; 10xobj: 6578
    end
    %     a = x(4); % um
    %     I0 = x(5);
    %     D = 2500;
    I = exp(-2*mus*(z+g.zs*3.03)) *(f^2)./ ( 4*(n*zp+f).^2 .* (1 + (pi*D^2./(4*lambda*(n*zp+f))).^2 .* (1-(n*zp+f)/(f)).^2) );
    %     I = exp(-2*mus*z) ./ ( 4*(zp+f).^2 .* (1 + (pi*D^2./(4*lambda*(zp+f))).^2 .* (1-(zp+f)/(f)).^2) );
    % %    Iold = exp(-2*mus*z) ./ ( 2*(n*zp+f).^2 .* (1 + (n*zp*pi*D^2./(4*lambda*f*(n*zp+f))).^2 ));
    
    if 0
        foo = (IoffSlope/Ioff)*z + 1;
        Io = sum(yy'.*I.*foo-Ioff*I.*foo.^2) / sum(I.^2.*foo.^2);
        I = (Io.*I + Ioff).*foo;
    else
        Io = sum(yy'.*I) / sum(I.^2);
        I = Io.*I;
    end
    g.Io = Io;
    
end