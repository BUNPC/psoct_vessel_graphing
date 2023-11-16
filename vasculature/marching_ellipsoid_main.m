%% Main script for performing marching ellipsoid
%{
This package was initially created by collaborators of David Boas. It has
been modified over the years. This main script is still a work in progress.

This script was run with a subset of a larger g. The code did not
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
    subdir = '\dist_corrected\volume\';
    sigdir = '\gsigma_1-3-5_gsize_5-13-21\';
    % Segmentation filename
    vol_name = 'ref_4ds_norm_inv_crop2.tif';
    % Graph filename
    graph_name = 'ref_4ds_norm_inv_crop2_segment_pmin_0.23_mask40_ds_mean_ds_graph.mat';
    graph_me_name = 'ref_4ds_norm_inv_crop2_segment_pmin_0.23_mask40_ds_mean_ds_marching_ellipse.mat';
%     graph_me_name = 'ref_4ds_norm_inv_crop2_segment_pmin_0.23_mask40_ds_mean_ds_mstep4_cutoff4';
    % View flag for marching ellipsoid (1 = view).
    viewflag = 1;
elseif isunix
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
    % Subject IDs
    subid = '/NC_6839/';
    subdir = '/dist_corrected/volume/';
    sigdir = '/gsigma_1-3-5_gsize_5-13-21/';
    % Segmentation filename
    vol_name = 'ref_4ds_norm_inv_crop2.tif';
    % Graph filename
    graph_name = 'ref_4ds_norm_inv_crop2_segment_pmin_0.23_mask40_ds_mean_ds_graph.mat';
    % View flag for marching ellipsoid (1 = view).
    viewflag = 0;
end

%% Load volume and g from Data
vox_dim = [12, 12, 15];

%%% Load volume
fullpath = fullfile(dpath, subid, subdir);
filename = strcat(fullpath, vol_name);
vol = mat2gray(TIFF2MAT(filename));
% volumeViewer(vol);

%%% Load graph
fullpath = fullfile(dpath, subid, subdir, sigdir);
filename = strcat(fullpath, graph_name);
g = load(filename);
g = g.im_re;
nodes = g.nodes;
edges = g.edges;
% Copy edges into standard format
s = edges(:,1); % source node
t = edges(:,2); % target node

% Create standard Matlab g
g_mat = graph(s, t);
figure;
p = plot(g_mat, 'XData', nodes(:,1), 'YData', nodes(:,2), 'ZData', nodes(:,3));
p.EdgeColor = 'red'; p.LineWidth = 1.5;
xlabel('x'); ylabel('y'); zlabel('z'); title('Before Marching Ellipsoid');
view(3);

%% Perform marching ellipsoid to connect disparate segments
flag_seeds = zeros(1, size(g.nodes,1));
graph2.nodes = [];
graph2.edges = [];
% cut off distance for neighboring node detection (voxels)
cutoff_dist = 4;
% Track number of segments that have been examined
i = 0;
ori = [];
% Increment (voxels) for each marching step
marching_step = 4;

%%% Iterate over each node in the graph
while sum(flag_seeds) < size(g.nodes,1)
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
    seed = round(g.nodes(idx(seed_idx),:));

    % Round the coordinates (this may be unecessary)
    cen_x = seed(1); cen_x = max(min(cen_x,size(vol,1)), 1);
    cen_y = seed(2); cen_y = max(min(cen_y,size(vol,2)), 1);
    cen_z = seed(3); cen_z = max(min(cen_z,size(vol,3)), 1);
    
    % fit an initial ellipsoid to that node, remember the sign of the
    % primary direction
    s = EllipseFit3DConstrained_jy(vol, cen_x, cen_y, cen_z, viewflag);
    a = [s.a1 s.a2 s.a3];
    [~,axis_idx] = max(a);
    vec = primary_dir(s,axis_idx);
        
    % remember the initial centroid
    ini_cen = s.mu;
    ini_vec = vec;
    ini_s = s;
    
    % flag for boundary detection
    boundary = 0;
    
    % record all seed nodes index
    if j == 0
        ori = [ori, size(graph2.nodes,1)+1];
    end

    %% Search positive direction
    disp('Start searching positive direction');
    % determine if it is an ellipsoid or not
    % if spheres are detected terminate the marching process
    while max(a)/min(a)>=1.2 && boundary==0 && vol(cen_x,cen_y,cen_z)>0.01
        j = j+1;
        % determine the maximum axis, travel along that axis
        [~,axis_idx] = max(a);
        vec_pre = vec;
        % first add the centroid into nodes list of new g
        graph2.nodes = [graph2.nodes; round(s.mu)'];
        nodes = unique(graph2.nodes,'rows','stable');
        
        % if reached local minimal
        if(size(nodes,1) ~= size(graph2.nodes,1))
            % repeat fitting on the same centroid & increase step size
            boundary = 1;
            graph2.nodes = nodes;
        else
            if size(graph2.nodes,1)>1 && j>1
                graph2.edges = [graph2.edges;...
                    [size(graph2.nodes,1)-1 size(graph2.nodes,1)] ];
            end
        end
        
        % if the controid reaches boundary, break the loop
        if cen_x>size(vol,1)-2 || cen_x<3 || cen_y>size(vol,2)-2 || cen_y<3 || cen_z>size(vol,3)-2 || cen_z<3
            boundary=1;
        end
        
        % determine the primary direction of fitted ellipsoid and move to
        % positive direction
        vec = primary_dir(s,axis_idx);
        angle = atan2(norm(cross(vec_pre,vec)), dot(vec_pre,vec));
        if(angle>pi/2)
            vec = -vec;
        end
        if abs(angle)>pi/6
            boundary = 1;
        end
        
        % flag all the nodes within this ellipsoid
        % find all adjacent nodes from vesselness filter g
        for num=1:length(idx)
            dist=sqrt(sum(g.nodes(idx(num),:) - [cen_x cen_y cen_z]).^2);
            if dist < cutoff_dist
                node=g.nodes(idx(num),:);
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
        cen_x=cen_new(1);cen_x=min(cen_x,size(vol,1));cen_x=max(cen_x,1);
        cen_y=cen_new(2);cen_y=min(cen_y,size(vol,2));cen_y=max(cen_y,1);
        cen_z=cen_new(3);cen_z=min(cen_z,size(vol,3));cen_z=max(cen_z,1);
        % update ellipsoid
        s = EllipseFit3DConstrained_jy(vol,cen_x,cen_y,cen_z, viewflag);
        % update the parameters
        a = [s.a1 s.a2 s.a3];
    end
    
    %% Search negative direction
    % Re-initialize boundary, vector
    boundary = 0;
    vec = -ini_vec;
    cen_new=round(ini_cen+marching_step.*vec);
    cen_x=cen_new(1);cen_x=min(cen_x,size(vol,1));cen_x=max(cen_x,1);
    cen_y=cen_new(2);cen_y=min(cen_y,size(vol,2));cen_y=max(cen_y,1);
    cen_z=cen_new(3);cen_z=min(cen_z,size(vol,3));cen_z=max(cen_z,1);
    
    % Run ellipsoid fitting
    s = EllipseFit3DConstrained_jy(vol,cen_x,cen_y,cen_z,viewflag);
    % Update local parameters
    a = [s.a1 s.a2 s.a3];
    neg=1;

    disp('Start searching negative direction');
    while max(a)/min(a)>=1.2 && boundary==0 && vol(cen_x,cen_y,cen_z)>0.01
        j=j+1;
        % determine the maximum axis, travel along that axis
        [~,axis_idx]=max(a);
        vec_pre=vec;
        % first add the centroid into nodes list of new g       
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
        if cen_x>size(vol,1)-2 || cen_x<3 || cen_y>size(vol,2)-2 || cen_y<3 || cen_z>size(vol,3)-2 || cen_z<3
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
            dist=sqrt(sum(g.nodes(idx(num),:) - [cen_x cen_y cen_z]).^2);
            if dist < cutoff_dist
                node=g.nodes(idx(num),:);
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
        cen_x=cen_new(1);cen_x=min(cen_x,size(vol,1));cen_x=max(cen_x,1);
        cen_y=cen_new(2);cen_y=min(cen_y,size(vol,2));cen_y=max(cen_y,1);
        cen_z=cen_new(3);cen_z=min(cen_z,size(vol,3));cen_z=max(cen_z,1);
        
        % update ellipsoid
        s = EllipseFit3DConstrained_jy(vol,cen_x,cen_y,cen_z,viewflag);
        % update the parameters
        a=[s.a1 s.a2 s.a3];     
    end
    disp(['Finish searching segment No.' num2str(i)]);
end
%%% Add results to struct
g = struct;
g.nodes = graph2.nodes;
g.edges = graph2.edges;

%% Remove single-point edges (self loops) from the graph
% Array for tracking single point edges
self_loop_idx = [];
for u = 1:size(g.edges,1)
    if g.edges(u,1) == g.edges(u,2)
        self_loop_idx = [self_loop_idx; u];
    end
end
% Remove self-loops from graph
g.edges(self_loop_idx,:) = [];

%% Remove floating nodes (without edges) from the graph
% Sometimes the graph from the output of the marching ellipsoid includes
% nodes that are not connected to any edges. In this case, it results in an
% error when converting it into a Matlab graph. The Matlab graph will
% expect fewer nodes than are provided. One solution is to determine if the
% last node in the list of nodes is a member of the edges array. If not,
% then remove this node.

%{
% Array of indices of all nodes
nidcs = 1:length(g.nodes);
% Find indices that are not in the edges array
ndel = find(~ismember(nidcs, g.edges) == 1);
% Remove nodes
g.nodes(ndel,:) = [];
%}

%% Save final graph
% Output directory
t = strcat('mstep',num2str(marching_step),'_cutoff',num2str(cutoff_dist));
graph_output = strcat(graph_name(1:end-9),t,'.mat');
fullpath = fullfile(dpath, subid, subdir, sigdir);
output_file = strcat(fullpath, graph_output);

% Save graph
save(output_file,'g');

function plot_graph(g, nodes, title_str)
figure;
p = plot(g, 'XData', nodes(:,1), 'YData', nodes(:,2), 'ZData', nodes(:,3));
p.EdgeColor = 'red'; p.LineWidth = 1.5;
xlabel('x'); ylabel('y'); zlabel('z'); title(title_str);
set(gca, 'FontSize', 20); grid on;
view(3);
end