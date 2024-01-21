

%%
function [I_nodes] = sk3D(sz,Graph,name,res,downsampled_flag,save_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [I_nodes] = sk3D([1100,1000,400],Data.Graph,'seg.sk.nii',[0.01,0.01,0.01],1,1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  sz     =  the size of output volume. The order is [ y x z ]
%
%            graph.nodes coordinates             has [ x y z ] order
%            nii size from mri_info              has [ x y z ] order
%
%            vesGraphValidate GUI image info     has [ y x z ] order
%            nii size after MRIread              has [ y x z ] order
%
%            Data.angio vol                      has [ z y x ] order
%  Graph  =  Graph struct that has node coordinate vector and edge vector 
%  name   =  filename to save as
%  res    =  resolution of the nifti image. 
%            i.e. [0.01, 0.01, 0.01] 0.01cm (10um) isotropics
%  downsampled_flag = 1=downsampled; 0=no-downsampled
%  save_flag = to save the skeleton or not
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% other function needed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  bresenham_line3d
%  MRIwrite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get a list of xyz coordinates of the nodes
if downsampled_flag==1  % 1 = graph has been downsampled
    edges = Graph.edges;
    nodes = round(Graph.nodes);
    
    pt1 = nodes(edges(:,1),:);
    pt2 = nodes(edges(:,2),:);
    
    % Since number of nodes were downsampled, bresenham_line3d func will
    % fill the grap between two adjacent nodes. 
    list = zeros(0,3);

    idx=find(sum((pt1-pt2).^2,2)<200^2); % find node pairs whose length does not exceed 200px (then exclude these artifact nodes)
    counts_artifact_edge = size(pt1,1)-size(idx,1);
    fprintf('find %d artifact edges\n',counts_artifact_edge);
    for i = idx.'
        list = [list; bresenham_line3d(pt1(i,:), pt2(i,:))];
    end
    list = [list;nodes];
else 
    edges = Graph.edges;
    nodes = round(Graph.nodes);
    
    pt1 = nodes(edges(:,1),:);
    pt2 = nodes(edges(:,2),:);
    
    idx=find(sum((pt1-pt2).^2,2)>3); % find connected node pairs that are not adjacent to each other

    list = zeros(0,3);
    for i = idx.'
        list = [list; bresenham_line3d(pt1(i,:), pt2(i,:))];
    end
    list = [list;nodes];
end

%%
I_nodes = get_sk(sz, list);
if save_flag==1
    save_mri(I_nodes,name,res,'uchar')
end
end

function I_nodes = get_sk(sz, list)

    list = list(:,[2,1,3]);                         % xyz -> yxz
    
    fprintf('max yxz = %i %i %i\n',max(list,[],1))  % could be the correct sz 
    fprintf('min yxz = %i %i %i\n',min(list,[],1))
    
    
    list(list==0)=1;                                % if coordinate is 0, change it to 1
    if any([list./sz]>1,"all")
        sprintf('max pos exceed vol size\n');
        list(list(:,1)>sz(1),1)=sz(1);
        list(list(:,2)>sz(2),2)=sz(2);
        list(list(:,3)>sz(3),3)=sz(3);
    end
    
    lnr = sub2ind(sz,list(:,1),list(:,2),list(:,3)); % y x z
    I_nodes = zeros(sz);
    I_nodes(lnr) = 1;
end

function save_mri(I, name, res, datatype)
    if ~exist("datatype","var")
        datatype = 'float';
    end

    disp(' - making hdr...');
    % Make Nifti and header
    colres = res(2); 
    rowres = res(1); 
    sliceres = res(3); 
    % mri.vol = I;
    mri.volres = [res(1) res(2) res(3)];
    mri.xsize = rowres;
    mri.ysize = colres;
    mri.zsize = sliceres;
    a = diag([-colres rowres sliceres 1]);
    mri.vox2ras0 = a;
    mri.volsize = size(I); if length(mri.volsize)==2;mri.volsize = [mri.volsize 1];end %% make 2D into 3D
    mri.vol = I;
    % mri.vol = flip(mri.vol,1);
    MRIwrite(mri,name,datatype);
    disp(' - done - ');
end