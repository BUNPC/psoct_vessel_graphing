function ShowLocalDataWithSE(Vol,s)
    [d1,d2,d3] = size(Vol);
    % b is the scaling factor of visualization
    b=10;%round(min([d1,d2,d3])/2)-10;

    Bl = round(s.mu - b);
    Bh = round(s.mu + b);
    v = zeros(2*b+1,2*b+1,2*b+1);
    for x=Bl(1):Bh(1)
        for y=Bl(2):Bh(2)
            for z=Bl(3):Bh(3)
                if (x>0)&&(y>0)&&(z>0)&&(x<=d2)&&(y<=d1)&&(z<=d3)
                    v(y-Bl(2)+1,x-Bl(1)+1,z-Bl(3)+1) = Vol(y,x,z);
                end
            end
        end
    end
    MIPshow(v);
    s.mu = s.mu - [Bl(1) Bl(2) Bl(3)]';
    renderEllipsoidonMIPs(1,s,2*b+1,2*b+1,2*b+1);
end