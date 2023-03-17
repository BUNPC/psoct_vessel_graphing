function Vec=primary_dir(s,axis_idx)
    % determine the primary direction of fitted ellipsoid
    XradiusVec = [[0 s.a1]' [0 0]' [0 0]']';
    YradiusVec = [[0 0]' [0 s.a2]' [0 0]']';
    ZradiusVec = [[0 0]' [0 0]' [0 s.a3]']';
    % calculated the direction based on quaternion s.q
    R = rotation_quat(s.q);
    XradVec = R*XradiusVec + repmat(s.mu,1,size(XradiusVec,2));
    YradVec = R*YradiusVec + repmat(s.mu,1,size(YradiusVec,2));
    ZradVec = R*ZradiusVec + repmat(s.mu,1,size(ZradiusVec,2));
    % these are vectors of ellipse radiuses translated to origin:
    Xvec = XradVec(:,2)-XradVec(:,1);
    Yvec = YradVec(:,2)-YradVec(:,1);
    Zvec = ZradVec(:,2)-ZradVec(:,1);
    if axis_idx==1
        Vec=Xvec./sqrt(sumsqr(Xvec));
    elseif axis_idx==2
        Vec=Yvec./sqrt(sumsqr(Yvec));
    elseif axis_idx==3
        Vec=Zvec./sqrt(sumsqr(Zvec));
    end
end