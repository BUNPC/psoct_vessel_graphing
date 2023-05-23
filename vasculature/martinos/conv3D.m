function [I_conv] = conv3D(I,kernel)

I_conv = zeros(size(I),'single');

for ii = 1:size(I,3)
    imageData=I(:,:,ii); % single(imread(filename(1).name,i));
    tmp=single(convn(imageData, ones(kernel,kernel)./kernel^2,'same')+0.0001);
    I_conv(:,:,ii)=single(imageData./tmp);
end
end