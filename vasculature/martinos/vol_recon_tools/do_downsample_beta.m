function I = do_downsample_beta(I,ds_xyz,method)
if nargin == 2; method = 'mean'; end

datatype = 'single';

sz = size(I);
sz_Id = floor(sz./ds_xyz);

dsd = find(ds_xyz~=1); %% find downsampling dimension
dimorder = {[1 2 3],[2 1 3],[3 1 2]};


for d = dsd
    I = permute(I,dimorder{d});

    ds = ds_xyz(d); % downsample factor
    
    % generate template
    sz(d) = floor(sz(d)/ds);         % size after downsample
    sz_perm = sz(dimorder{d});       % size after downsample and permutation
    tmp = zeros(sz_perm,datatype);   % template vol after ds and perm 

    % perform downsample
    matrix = reshape(1:sz_perm(1)*ds,ds,[]);
    for i = 1:sz_perm(1)             % each plane in template vol
        switch method
            case 'mean'
                tmp(i,:,:) = mean(I(matrix(:,i),:,:),1);
            case 'max'
                tmp(i,:,:) = max(I(matrix(:,i),:,:),[],1);
            case 'min'
                tmp(i,:,:) = mminean(I(matrix(:,i),:,:),[],1);
        end
    end


    I = ipermute(I,dimorder{d});     % inverse permutation
end
    

%%

% tmp = zeros(sz_Id,'half');
% 
% Id = zeros(sz_Id,'single');
% 
% d3_mtx = reshape(1:sz_Id(3)*ds_xyz(3),ds_xyz(3),[]);
% 
% for d3 = 1:sz_Id(3)
%     tic_mean = tic;
%     fprintf('%i ',d3);
%     Id(:,:,d3) = mean(single(Ma(:,:,d3_mtx(1,d3):d3_mtx(end,d3))),3);
%     toc_mean=toc(tic_mean);
%     fprintf('%f seconds \n',toc_mean);
% end
end