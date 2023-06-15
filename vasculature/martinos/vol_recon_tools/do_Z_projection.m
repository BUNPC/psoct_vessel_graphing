function IP = do_Z_projection(ds,vol,method)
sz = size(vol);
fprintf('project range = %i\n', ds)
fprintf('volume size = %i %i %i\n', sz)

tic
for i = 1:ds
    I = do_downsample(vol(:,:,i:end-1-ds+i),[1,1,ds],method);
    if i==1; IP = I;
    else;    IP = cat(2,IP,I); 
    end
end

IP = reshape(IP,sz(1),sz(2),[]);
switch method
    case 'min'
        I_last = min(vol(:,:,end-ds+1:end),[],3);
    case 'max'
        I_last = max(vol(:,:,end-ds+1:end),[],3);
end

sz2 = size(IP);

IP = cat(3,IP, repmat(I_last,[1,1,sz(3)-sz2(3)]));
toc
end
