%% Main script for performing marching ellipsoid
%{
This package was initially created by collaborators of David Boas. It has
been modified over the years. This main script is still a work in progress.

This script was run with a subset of a larger graph. The code did not
compute the graph2.edges. Need to debug this.
%}
clear; clc; close all;
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
    subdir = '\dist_corrected\volume\gsigma_1-2-3-4-5_gsize_5--9-13-17-21\';
    % Segmentation filename
    seg_name = 'ref_4ds_norm_inv_segment_pmin_0.23_mask40_crop';
    % filename extension
    ext = '.tif';
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
vox_dim = [12, 12, 15];

%%% Convert segment to graph
fullpath = fullfile(dpath, subid, subdir);
filename = strcat(fullpath, strcat(seg_name, ext));
seg = TIFF2MAT(filename);
graph = seg_to_graph(seg, vox_dim);

%%% Load segmentation and graph
%{
% Define entire filepath 
fullpath = fullfile(dpath, subid, subdir);
filename = strcat(fullpath, strcat(fname, ext));

% Load Data
Data = load(filename, 'Data');

% Extract angio (segmentation (uint8)) from Data
seg = Data.Data.angio;
graph = Data.Data.Graph;

% Reorder the angio
seg = permute(seg, [3,2,1]);
% Normalize segmentation
seg = (seg-min(seg(:)))./(max(seg(:))-min(seg(:)));
%}
%% Generate seed points coordinates from intensity threshold
% This section was commented out. It appears to find voxels within the
% probability map that are above a threshold. The threshold is an array
% though, so the line "k = find(seg(:) > thresh)" crashed because there is
% not enough memory. I believe this section was just for testing and not
% implementation.

%{
% Initialize threshold array
min_th=0.5;
max_th=0.55;
th_step=(max_th-min_th)/64;
thresh = min_th : th_step : max_th;
idx1=[]; idx2=[]; idx3=[];

% Find voxels greater than threshold
V_bi = zeros(size(seg));
for z=1:size(seg,3)
    rl_depth = mod(z-1,65)+1;
    tmp = squeeze(seg(:,:,z));
    tmp2 = zeros(size(tmp));
    tmp = (tmp-min(tmp(:)))./(max(tmp(:))-min(tmp(:)));
    if(z==1)
        k = find(tmp(:)>0.3);
    else
        k = find(tmp(:)>thresh(rl_depth));
    end
    tmp2(k)=1;
    % m(z)=length(k);
    [tidx1, tidx2] = ind2sub(size(tmp),k);
    V_bi(:,:,z) = tmp2;
    idx1 = [idx1; tidx1];
    idx2 = [idx2; tidx2];
    idx3 = [idx3; ones(length(tidx1),1).*z];
end

% Output the seed points
% MAT2TIFF(V_bi,'th_bi.tif');
k = find(seg(:) > thresh);
[idx1,idx2,idx3] = ind2sub(size(seg),k);
graph.nodes = [idx1,idx2,idx3];
%}
%% perform marching ellipsoid to generate a new graph
flag_seeds = zeros(1, size(graph.nodes,1));
graph2.nodes = [];
graph2.edges = [];
% cut off distance for neighboring node detection (voxels)
cutoff_dist = 15;
% Track number of segments that have been examined
i = 0;
ori = [];

% for visualization purpose
% h=figure;
% x0=100;
% y0=100;
% width=1000;
% height=800;
% set(gcf,'position',[x0,y0,width,height]);
% axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'marchinging3.gif';

%%% Iterate over each node in the graph
while sum(flag_seeds) < size(graph.nodes,1)
    %%% pick a random node as seed point
    % Iterate number of segments tested
    i = i+1;
    % number of marching ellipsoids
    j = 0;
    % Find indices where flag_seeds equals zero
    [idx,~] = find(flag_seeds(:)==0);
    % Choose random integer within [1, length(idx)]
    seed_idx = randi(length(idx),1);
    % Set seed flag equal to 1 for these random integer indices
    flag_seeds(idx(seed_idx)) = 1;
    
    % seed = node coordinates of the seed_idx
    seed = round(graph.nodes(idx(seed_idx),:));

    % Round the coordinates (this may be unecessary)
    cen_x = seed(1); cen_x = max(min(cen_x,size(seg,1)), 1);
    cen_y = seed(2); cen_y = max(min(cen_y,size(seg,2)), 1);
    cen_z = seed(3); cen_z = max(min(cen_z,size(seg,3)), 1);
    
    % marching step size, need to be tuned in the future
    marching_step = 8;
    
    % fit an initial ellipsoid to that node, remember the sign of the
    % primary direction
    s = EllipseFit3DConstrained_jy(seg,cen_x,cen_y,cen_z,0);
    a = [s.a1 s.a2 s.a3];    
    [~,axis_idx] = max(a);
    vec=primary_dir(s,axis_idx);
    
%     % show the data & save to gif
%     ShowLocalDataWithSE(seg,s);
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
        ori = [ori, size(graph2.nodes,1)+1];
    end
    
    disp('Start searching positive direction');
    
    % determine if it is an ellipsoid or not
    % if spheres are detected terminate the marching process
    while max(a)/min(a)>=1.2 && boundary==0 && seg(cen_x,cen_y,cen_z)>0.01
        j=j+1;
        % determine the maximum axis, travel along that axis
        [~,axis_idx]=max(a);
        vec_pre=vec;
        % first add the centroid into nodes list of new graph
        graph2.nodes=[graph2.nodes; round(s.mu)'];
        nodes = unique(graph2.nodes,'rows','stable');
        % if reached local minimal (repeating fitting on the same
        % conetroid), increase step size
        if(size(nodes,1)~=size(graph2.nodes,1))
            boundary=1;
            graph2.nodes=nodes;
        else
            if size(graph2.nodes,1)>1 && j>1
                graph2.edges=[graph2.edges;[size(graph2.nodes,1)-1 size(graph2.nodes,1)]];
            end
        end
        
        
        % if the controid reaches boundary, break the loop
        if cen_x>size(seg,1)-2 || cen_x<3 || cen_y>size(seg,2)-2 || cen_y<3 || cen_z>size(seg,3)-2 || cen_z<3
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
            dist=sqrt(sum(graph.nodes(idx(num),:) - [cen_x cen_y cen_z]).^2);
            if dist < cutoff_dist
                node=graph.nodes(idx(num),:);
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
        cen_x=cen_new(1);cen_x=min(cen_x,size(seg,1));cen_x=max(cen_x,1);
        cen_y=cen_new(2);cen_y=min(cen_y,size(seg,2));cen_y=max(cen_y,1);
        cen_z=cen_new(3);cen_z=min(cen_z,size(seg,3));cen_z=max(cen_z,1);
        
        % update ellipsoid
        s = EllipseFit3DConstrained_jy(seg,cen_x,cen_y,cen_z,0);
%         ShowLocalDataWithSE(seg,s);
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
    cen_x=cen_new(1);cen_x=min(cen_x,size(seg,1));cen_x=max(cen_x,1);
    cen_y=cen_new(2);cen_y=min(cen_y,size(seg,2));cen_y=max(cen_y,1);
    cen_z=cen_new(3);cen_z=min(cen_z,size(seg,3));cen_z=max(cen_z,1);
    disp('Start searching negative direction');
    
    s = EllipseFit3DConstrained_jy(seg,cen_x,cen_y,cen_z,0);
    a=[s.a1 s.a2 s.a3];
    neg=1;
    
    while max(a)/min(a)>=1.2 && boundary==0 && seg(cen_x,cen_y,cen_z)>0.01
        j=j+1;
        % determine the maximum axis, travel along that axis
        [~,axis_idx]=max(a);
        vec_pre=vec;
        % first add the centroid into nodes list of new graph       
        graph2.nodes=[graph2.nodes; round(s.mu)'];
        nodes = unique(graph2.nodes,'rows','stable');
        % if reached local minimal (repeating fitting on the same
        % conetroid), increase step size
        if(size(nodes,1)~=size(graph2.nodes,1))
            boundary=1;
            graph2.nodes=nodes;
        else
            if(neg==1)
                if j>1
                    graph2.edges=[graph2.edges;[ori(i) size(graph2.nodes,1)]];
                end
                neg=0;
            else
                if j>1
                    graph2.edges=[graph2.edges;[size(graph2.nodes,1)-1 size(graph2.nodes,1)]];
                end
            end
        end
        
        
        % if the controid reaches boundary, break the loop
        if cen_x>size(seg,1)-2 || cen_x<3 || cen_y>size(seg,2)-2 || cen_y<3 || cen_z>size(seg,3)-2 || cen_z<3
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
            dist=sqrt(sum(graph.nodes(idx(num),:) - [cen_x cen_y cen_z]).^2);
            if dist < cutoff_dist
                node=graph.nodes(idx(num),:);
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
        cen_x=cen_new(1);cen_x=min(cen_x,size(seg,1));cen_x=max(cen_x,1);
        cen_y=cen_new(2);cen_y=min(cen_y,size(seg,2));cen_y=max(cen_y,1);
        cen_z=cen_new(3);cen_z=min(cen_z,size(seg,3));cen_z=max(cen_z,1);
        
        % update ellipsoid
        s = EllipseFit3DConstrained_jy(seg,cen_x,cen_y,cen_z,0);
%         ShowLocalDataWithSE(seg,s);
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

graph = [];
graph.nodes = graph2.nodes;
graph.edges = graph2.edges;

%% remove unconnected components in the graph
sameEdgeIdx = [];
for u = 1:size(graph.edges,1)
    if graph.edges(u,1) == graph.edges(u,2)
        sameEdgeIdx = [sameEdgeIdx; u];
    end
end
graph.edges(sameEdgeIdx,:) = [];

%% save the final graph
save('Graph_ellipsoid.mat','graph');

