function [Xvec,Yvec,Zvec]=dir_vec(s)
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
end