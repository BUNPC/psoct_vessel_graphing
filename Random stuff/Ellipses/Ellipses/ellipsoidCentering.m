function [nodePos] = ellipsoidCentering(I,nodePos,fillFactor)


nodeX=nodePos(:,1);
nodeY=nodePos(:,2);
nodeZ=nodePos(:,3);

movedNodesFlag=zeros(size(nodePos,1),1); % 1 = moved; 0 = not moved

hwait=waitbar(0,'Generating ellipsoids and centering nodes...');



for node=1:size(nodePos,1)
    if movedNodesFlag(node)==0
        waitbar(node/size(nodePos,1),hwait);

        s = EllipseFit3DConstrained_dab( I,nodeX(node),nodeY(node),nodeZ(node),0); % the '0' turns off the 'viewFlag'

        % how to find points close to ellipsoid center, inside a cylinder defined
        % by ellipsoid radius and height equal to fill factor length

        % fill factor = f (it should be hard coded somewhere in imView3d - length
        % between nodes; jmoose: I chose 'fillFactor' from the fill function

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

        % these are just vectors of ellipse radiuses translated to origin:
        Xvec = XradVec(:,2)-XradVec(:,1);
        Yvec = YradVec(:,2)-YradVec(:,1);
        Zvec = ZradVec(:,2)-ZradVec(:,1);
        
        Xvec=Xvec/norm(Xvec);
        Yvec=Yvec/norm(Yvec);
        Zvec=Zvec/norm(Zvec);
        % They can be generated faster than in above code, but lets not think about this for now...

        % for any node at point (X,Y,Z) in vessel stack, we can estimate if it is
        % inside cylinder around ellipse's center...

        % loop through all nodes in graph...
        for node2=1:size(nodePos,1)
            
            if movedNodesFlag(node2)==0 % should hopefully get faster exponentially each time through first 'for' loop
        
                distFromEllipCenter = nodePos(node2,:)' - s.mu; % example signle node (x, y, z) coordinates are [2, 4, 6]. Needed to subtract s.mu to get node in respect to ellipse center

                if (abs(distFromEllipCenter'*Yvec)<norm(fillFactor)) && ( abs( distFromEllipCenter'*Xvec + distFromEllipCenter'*Zvec )< max(s.a1, s.a3))
                   % looks at scalar vector products to decide if point is inside cylinder...
                   % point is inside cylinder, move the node and set its flag
                   nodePos(node2,:)=s.mu';
                   movedNodesFlag(node2) = 1;
                end;
            end
        end
    end
end


close(hwait);
