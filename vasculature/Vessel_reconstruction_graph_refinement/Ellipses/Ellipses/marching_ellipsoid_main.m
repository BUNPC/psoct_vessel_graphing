%% Main script for performing marching ellipsoid
%{
This package was initially created by collaborators of David Boas. It has
been modified over the years. This main script is still a work in progress.
%}

%% Add top-level directory of code repository to path
% This allows Matlab to find the functions in the project folders

% Start in current directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Truncate path to reach top-level directory (psoct_vessel_graphing)
topdir = mydir(1:idcs(end));
addpath(genpath(topdir));

%% Initialize data path for linux or personal machine (debugging)

%%% Local machine
if ispc
    dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures\test_data\Ann_Mckee_samples_10T\';
    % Subject IDs
    subid = 'NC_6839';
    subdir = '\dist_corrected\volume\gsigma_7--9-11_gsize_29-37-45\';
    % Filename to parse (this is test data)
    fname = 'ref_4ds_norm_inv_segment_pmin_0.23_mask_40_graph_data';
    % filename extension
    ext = '.mat';
%%% Computing cluster (SCC)
elseif isunix
    % Path to top-level directory
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_10T/';
    % Subfolder containing data
    subdir = '/dist_corrected/volume/';
    % Filename to parse (this will be the same for each subject)
    fname = 'ref_4ds_norm_inv';
    % filename extension
    ext = '.tif';
    % Complete subject ID list for Ann_Mckee_samples_10T
    subid = {'AD_10382', 'AD_20832', 'AD_20969', 'AD_21354', 'AD_21424',...
             'CTE_6489', 'CTE_6912', 'CTE_7019', 'CTE_8572', 'CTE_7126',...
             'NC_21499', 'NC_6047', 'NC_6839', 'NC_6974', 'NC_7597',...
             'NC_8095', 'NC_8653'};
        
    %%% Create cell array of subject ID and sigma for job array on the SCC 
    nrow = length(subid)*size(sigmas,2);
    nsigma = size(sigmas,2);
    sub_sigma = cell(length(subid).*size(sigmas,2), 2);
    idx = 1;
    % Fill sub_sigma cell array with each sigma array for each subject
    for i = 1:3:nrow
        sub_sigma{i,1} = subid{idx};
        sub_sigma{(i+1),1} = subid{idx};
        sub_sigma{(i+2),1} = subid{idx};
        idx = idx + 1;
        for j = 1:nsigma
            sub_sigma{(i+j-1),2} = sigmas(j,:);
        end
    end
    
    %%% Reassign subid and sigma based on job array counter
    % Retrieve SGE_TASK_ID from system (job array index)
    batch_idx = getenv('SGE_TASK_ID');
    
    % If this is a job array, then batch_idx will not be empty.
    if ~isempty(batch_idx)
        % Convert from ASCII to double
        batch_idx = str2double(batch_idx);
        % Retrieve corresponding row from sub_sigma
        [subid, gsigma] = sub_sigma{batch_idx, :};
    % Otherwise, set the Gaussian sigma manually
    else
        subid = 'NC_6839';
        gsigma = [7, 9, 11];
    end
end

%% Load segmentation volume and graph from Data
% Define entire filepath 
fullpath = fullfile(dpath, subid, subdir);
filename = strcat(fullpath, strcat(fname, ext));

% Load Data
Data = load(filename, 'Data');

% Extract angio (segmentation (uint8)) from Data
seg = Data.Data.angio;
graph = Data.Data.Graph;

% Normalize segmentation
seg = (seg-min(seg(:)))./(max(seg(:))-min(seg(:)));

%% generate coordinates of seed points based on intensity threshold
%{
% apply threshold for each xy slice
min_th=0.5;
max_th=0.55;
th_step=(max_th-min_th)/64;
thresh=min_th:th_step:max_th;
idx1=[];
idx2=[];
idx3=[];
V_bi=zeros(size(V));
for z=1:size(V,3)
    rl_depth=mod(z-1,65)+1;
    tmp=squeeze(V(:,:,z));
    tmp2=zeros(size(tmp));
    tmp=(tmp-min(tmp(:)))./(max(tmp(:))-min(tmp(:)));
    if(z==1)
        k=find(tmp(:)>0.3);
    else
        k=find(tmp(:)>thresh(rl_depth));
    end
    tmp2(k)=1;
    % m(z)=length(k);
    [tidx1,tidx2]=ind2sub(size(tmp),k);
    V_bi(:,:,z)=tmp2;
    idx1=[idx1; tidx1];
    idx2=[idx2; tidx2];
    idx3=[idx3; ones(length(tidx1),1).*z];
end
MAT2TIFF(V_bi,'th_bi.tif');
k=find(V(:)>thresh);
[idx1,idx2,idx3]=ind2sub(size(V),k);
Graph.nodes=[idx1,idx2,idx3];
%}
%% perform marching ellipsoid to generate a new graph
flag_seeds=zeros(1,size(Graph.nodes,1));
Graph_new.nodes=[];
Graph_new.edges=[];
cutoff_dist=15; % cut off distance for neighboring node detection

% for visualization purpose
% h=figure;
% x0=100;
% y0=100;
% width=1000;
% height=800;
% set(gcf,'position',[x0,y0,width,height]);
% axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'marchinging3.gif';

i=0;    % number of segments
ori=[];

while sum(flag_seeds)<size(Graph.nodes,1)
    
    % pick a random node as seed point
    i=i+1;  % keep tracking of number of segments
    j=0;    % number of marching ellipsoids
    [idx,~]=find(flag_seeds(:)==0);
    seed_idx=randi(length(idx),1);
    % set seed flag to be 1
    flag_seeds(idx(seed_idx))=1;
    seed=round(Graph.nodes(idx(seed_idx),:));
    cen_x=seed(1);cen_x=min(cen_x,size(V,1));cen_x=max(cen_x,1);
    cen_y=seed(2);cen_y=min(cen_y,size(V,2));cen_y=max(cen_y,1);
    cen_z=seed(3);cen_z=min(cen_z,size(V,3));cen_z=max(cen_z,1);
    
    marching_step=8;    % marching step size, need to be tuned in the future
    
    % fit an initial ellipsoid to that node, remember the sign of the
    % primary direction
    s = EllipseFit3DConstrained_jy(V,cen_x,cen_y,cen_z,0);
    a=[s.a1 s.a2 s.a3];    
    [~,axis_idx]=max(a);
    vec=primary_dir(s,axis_idx);
    
%     % show the data & save to gif
%     ShowLocalDataWithSE(V,s);
%     title(['Segment No.: ',num2str(i), ' Marching point No.: ',num2str(j)]);
%     frame = getframe(h); 
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256); 
%     % Write to the GIF File 
%     if i==1
%         imwrite(imind,cm,filename,'gif', 'Loopcount',Inf,'DelayTime',1/2); 
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/2);
%     end
        
    % remember the initial centroid
    ini_cen=s.mu;
    ini_vec=vec;
    ini_s=s;
    
    % flag for boundary detection
    boundary = 0;
    
    % record all seed nodes index
    if j==0
        ori=[ori size(Graph_new.nodes,1)+1];
    end
    
    disp('Start searching positive direction');
    
    % determine if it is an ellipsoid or not
    % if spheres are detected terminate the marching process
    while max(a)/min(a)>=1.2 && boundary==0 && V(cen_x,cen_y,cen_z)>0.01
        j=j+1;
        % determine the maximum axis, travel along that axis
        [~,axis_idx]=max(a);
        vec_pre=vec;
        % first add the centroid into nodes list of new graph
        Graph_new.nodes=[Graph_new.nodes; round(s.mu)'];
        nodes = unique(Graph_new.nodes,'rows','stable');
        % if reached local minimal (repeating fitting on the same
        % conetroid), increase step size
        if(size(nodes,1)~=size(Graph_new.nodes,1))
            boundary=1;
            Graph_new.nodes=nodes;
        else
            if size(Graph_new.nodes,1)>1 && j>1
                Graph_new.edges=[Graph_new.edges;[size(Graph_new.nodes,1)-1 size(Graph_new.nodes,1)]];
            end
        end
        
        
        % if the controid reaches boundary, break the loop
        if cen_x>size(V,1)-2 || cen_x<3 || cen_y>size(V,2)-2 || cen_y<3 || cen_z>size(V,3)-2 || cen_z<3
            boundary=1;
        end
        
        % determine the primary direction of fitted ellipsoid and move to
        % positive direction
        vec=primary_dir(s,axis_idx);
        angle = atan2(norm(cross(vec_pre,vec)), dot(vec_pre,vec));
        if(angle>pi/2)
            vec=-vec;
        end
        if abs(angle)>pi/6
            boundary=1;
        end
        
        % flag all the nodes within this ellipsoid
        % find all adjacent nodes from vesselness filter graph
        for num=1:length(idx)
            dist=sqrt(sum(Graph.nodes(idx(num),:) - [cen_x cen_y cen_z]).^2);
            if dist < cutoff_dist
                node=Graph.nodes(idx(num),:);
                [Xvec,Yvec,Zvec]=dir_vec(s);
                if abs((node- [cen_x cen_y cen_z])*Yvec)<a(2) && abs((node- [cen_x cen_y cen_z])*Xvec)<a(1) && abs((node- [cen_x cen_y cen_z])*Zvec)<a(3)
                    % looks at scalar vector products to decide if point is inside cylinder...
                    % point is inside cylinder, move the node and set its flag
                    flag_seeds(idx(num))=1;
                    % disp('node captured!');
                end
            end
        end
        
        % move to next centroid along this direction
        cen_new=round(s.mu+marching_step.*vec);
        cen_x=cen_new(1);cen_x=min(cen_x,size(V,1));cen_x=max(cen_x,1);
        cen_y=cen_new(2);cen_y=min(cen_y,size(V,2));cen_y=max(cen_y,1);
        cen_z=cen_new(3);cen_z=min(cen_z,size(V,3));cen_z=max(cen_z,1);
        
        % update ellipsoid
        s = EllipseFit3DConstrained_jy(V,cen_x,cen_y,cen_z,0);
%         ShowLocalDataWithSE(V,s);
%         title(['Segment No.: ',num2str(i), ' Marching point No.: ',num2str(j)]);
%         frame = getframe(h); 
%         im = frame2im(frame); 
%         [imind,cm] = rgb2ind(im,256); 
%         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/2);
        % update the parameters
        a=[s.a1 s.a2 s.a3];     
    end
    
    % go back to the starting point and march to negative direction
    boundary=0;
    vec=-ini_vec;
    marching_step=8;
    cen_new=round(ini_cen+marching_step.*vec);
    cen_x=cen_new(1);cen_x=min(cen_x,size(V,1));cen_x=max(cen_x,1);
    cen_y=cen_new(2);cen_y=min(cen_y,size(V,2));cen_y=max(cen_y,1);
    cen_z=cen_new(3);cen_z=min(cen_z,size(V,3));cen_z=max(cen_z,1);
    disp('Start searching negative direction');
    
    s = EllipseFit3DConstrained_jy(V,cen_x,cen_y,cen_z,0);
    a=[s.a1 s.a2 s.a3];
    neg=1;
    
    while max(a)/min(a)>=1.2 && boundary==0 && V(cen_x,cen_y,cen_z)>0.01
        j=j+1;
        % determine the maximum axis, travel along that axis
        [~,axis_idx]=max(a);
        vec_pre=vec;
        % first add the centroid into nodes list of new graph       
        Graph_new.nodes=[Graph_new.nodes; round(s.mu)'];
        nodes = unique(Graph_new.nodes,'rows','stable');
        % if reached local minimal (repeating fitting on the same
        % conetroid), increase step size
        if(size(nodes,1)~=size(Graph_new.nodes,1))
            boundary=1;
            Graph_new.nodes=nodes;
        else
            if(neg==1)
                if j>1
                    Graph_new.edges=[Graph_new.edges;[ori(i) size(Graph_new.nodes,1)]];
                end
                neg=0;
            else
                if j>1
                    Graph_new.edges=[Graph_new.edges;[size(Graph_new.nodes,1)-1 size(Graph_new.nodes,1)]];
                end
            end
        end
        
        
        % if the controid reaches boundary, break the loop
        if cen_x>size(V,1)-2 || cen_x<3 || cen_y>size(V,2)-2 || cen_y<3 || cen_z>size(V,3)-2 || cen_z<3
            boundary=1;
        end
        
        % determine the primary direction of fitted ellipsoid and move to
        % positive direction
        vec=primary_dir(s,axis_idx);
        angle = atan2(norm(cross(vec_pre,vec)), dot(vec_pre,vec));
        if(angle>pi/2)
            vec=-vec;
        end
        if abs(angle)>pi/6
            boundary=1;
        end
        
        % flag all the nodes within this ellipsoid
        for num=1:length(idx)
            dist=sqrt(sum(Graph.nodes(idx(num),:) - [cen_x cen_y cen_z]).^2);
            if dist < cutoff_dist
                node=Graph.nodes(idx(num),:);
                [Xvec,Yvec,Zvec]=dir_vec(s);
                if abs((node- [cen_x cen_y cen_z])*Yvec)<a(2) && abs((node- [cen_x cen_y cen_z])*Xvec)<a(1) && abs((node- [cen_x cen_y cen_z])*Zvec)<a(3)
                    % looks at scalar vector products to decide if point is inside cylinder...
                    % point is inside cylinder, move the node and set its flag
                    flag_seeds(idx(num))=1;
                    % disp('node captured!');
                end
            end
        end
        
        % move to next centroid along this direction
        cen_new=round(s.mu+marching_step.*vec);
        cen_x=cen_new(1);cen_x=min(cen_x,size(V,1));cen_x=max(cen_x,1);
        cen_y=cen_new(2);cen_y=min(cen_y,size(V,2));cen_y=max(cen_y,1);
        cen_z=cen_new(3);cen_z=min(cen_z,size(V,3));cen_z=max(cen_z,1);
        
        % update ellipsoid
        s = EllipseFit3DConstrained_jy(V,cen_x,cen_y,cen_z,0);
%         ShowLocalDataWithSE(V,s);
%         title(['Segment No.: ',num2str(i), ' Marching point No.: ',num2str(j)]);
%         frame = getframe(h); 
%         im = frame2im(frame); 
%         [imind,cm] = rgb2ind(im,256); 
%         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/2);
        % update the parameters
        a=[s.a1 s.a2 s.a3];     
    end
    disp(['Finish searching segment No.' num2str(i)]);
end

Graph=[];
Graph.nodes=Graph_new.nodes;
Graph.edges=Graph_new.edges;

%% remove unconnected components in the graph
sameEdgeIdx = [];
for u = 1:size(Graph.edges,1)
    if Graph.edges(u,1) == Graph.edges(u,2)
        sameEdgeIdx = [sameEdgeIdx; u];
    end
end
Graph.edges(sameEdgeIdx,:) = [];

%% save the final graph
save('Graph_ellipsoid.mat','Graph');

