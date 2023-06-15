function [Id] = do_downsample(I, ds_xyz, method)
if nargin == 2; method = 'mean'; end

ds_x=ds_xyz(1);
ds_y=ds_xyz(2);
ds_z=ds_xyz(3);

sz = size(I);

m=1:ds_y:sz(1)-ds_y+1; %dim1
n=1:ds_x:sz(2)-ds_x+1; %dim2
p=1:ds_z:sz(3)-ds_z+1; %dim3

dim1              = repmat(repmat(m.',[1, length(n)]),[1,1,length(p)]);
dim2              = repmat(repmat(n,[length(m), 1]),[1,1,length(p)]);
dim3              = reshape(ones(length(m)*length(n),1)*p, [length(m),length(n),length(p)]);

lnr_anchor          = sub2ind(sz, dim1, dim2, dim3);
sz_lnr_anchor       = size(lnr_anchor);
                    lnr_anchor = lnr_anchor(:); 

d1                = reshape([1:ds_y].' .* ones(1,ds_x*ds_z),[ds_y, ds_x, ds_z]); %[d1 d2 d3]
d1                = permute(d1, [1,2,3]);

d2                = reshape([1:ds_x].' .* ones(1,ds_y*ds_z),[ds_x, ds_y, ds_z]); %[d2 d1 d3]
d2                = permute(d2, [2,1,3]);

d3                = reshape([1:ds_z].' .* ones(1,ds_y*ds_x),[ds_z, ds_y, ds_x]); %[d3 d1 d2]
d3                = permute(d3, [2,3,1]);

lnr_box_mold        = sub2ind(sz,d1,d2,d3);
                    lnr_box_mold = lnr_box_mold(:).'; 

lnr_box             = lnr_anchor - 1 + lnr_box_mold;
sz_lnr_box          = size(lnr_box);

switch method
    case 'mean'
        Id = reshape(mean(reshape(I(lnr_box),sz_lnr_box),2),sz_lnr_anchor);
    case 'min'
        Id = reshape(min(reshape(I(lnr_box),sz_lnr_box),[],2),sz_lnr_anchor);
    case 'max'
        Id = reshape(max(reshape(I(lnr_box),sz_lnr_box),[],2),sz_lnr_anchor);
end
end