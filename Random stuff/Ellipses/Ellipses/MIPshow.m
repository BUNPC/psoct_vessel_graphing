function MIPshow(Vol)
    [d1,d2,d3] = size(Vol);
    %Myx = feval(F,Vol,[],3);
    Myx = max(Vol,[],3);
    %Mxz = feval(F,Vol,[],1);
    Mxz = max(Vol,[],1);
    Mxz = reshape(Mxz,d2,d3);
    %Myz = feval(F,Vol,[],2);
    Myz = max(Vol,[],2);
    Myz = reshape(Myz,d1,d3);
    % M = [imrotate(mat2gray(Myz),90) zeros(d3,d3); mat2gray(Mxy) (mat2gray(Mxz))];
    M = [Myx, Myz ; Mxz', zeros(d3,d3) ];
    imagesc(M, [0, 1]);
    colormap(gray);
    ylabel('Y (Dim 1)                                  Z (Dim 3)');
    xlabel('X (Dim 2)                                  Y (Dim 1)');
    % set(gca,'Ydir','normal');
end