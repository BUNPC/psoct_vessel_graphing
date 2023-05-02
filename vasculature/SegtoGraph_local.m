%% Convert the vessel segments to a graph (nodes + vertices)
%{
This script takes the output of the script ~\vasculature\vesSegment\vesSegment
and then creates a graph.
%}

%% Add top-level directory of code repository to path
% Print current working directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Remove the two sub folders to reach parent
% (psoct_human_brain\vasculature\vesSegment)
newdir = mydir(1:idcs(end));
addpath(genpath(newdir));

%% Hardcoded file for debugging:
dpath = 'C:\Users\mack\Documents\BU\Boas_Lab\psoct_human_brain_resources\test_data\Hui_Frangi_dataset\200218depthnorm';    % contained in subfolder
filename = '\volume_ori_inv_cropped_sigma2.mat';

%% Create the local vessel_mask variable
[~,~,ext] = fileparts(filename);
if strcmp(ext,'.mat')
    temp = load([dpath filename]);
    fn = fieldnames(temp);
    angio = temp.(fn{1});
elseif  strcmp(ext,'.tiff') || strcmp(ext,'.tif')
    info = imfinfo([dpath filename]);
    for u = 1:length(info)
        if u == 1
            temp = imread([dpath filename],1);
            angio = zeros([size(temp) length(info)]);
            angio(:,:,u) = temp;
        else
            angio(:,:,u) = imread([dpath filename],u);
        end
    end
end

% Convert all nonzero values into a logical 1.
vessel_mask = logical(angio);

% remove islands from raw segments
CC = bwconncomp(vessel_mask);
for uuu = 1:length(CC.PixelIdxList)
    if length(CC.PixelIdxList{uuu}) < 100    % 30 for default
        vessel_mask(CC.PixelIdxList{uuu}) = 0;
    end
end

%% Create graph from binary vessel mask
% Reduce 3-D binary volume to a curve skeleton
vessel_skl = bwskel(vessel_mask,'MinBranchLength',1);
% Compute graph from skeleton
vessel_graph = fun_skeleton_to_graph(vessel_skl);
% Compute Euclidean distance transform of inverse of binary vessel mask
vessel_mask_dt = bwdist(~vessel_mask);

%% Count the number of edges

% Size of the angiogram. It will help to convert indeces to subscripts
angio_size = size(vessel_mask);

% find length of nodes to allocate size
nodes_ind_count = length(vessel_graph.node.cc_ind)+length(vessel_graph.link.pos_ind);
nodes_ind = zeros(1,nodes_ind_count);

% find length of edges to allocate size
edges_ind_count = 0;
% link_cc_ind = zeros(size(vessel_graph.link.cc_ind));
for u = 1:length(vessel_graph.node.connected_link_label)
    edges_ind_count = edges_ind_count+length(vessel_graph.node.connected_link_label{u});
end
for  u = 1:length(vessel_graph.link.cc_ind)
    edges_ind_count = edges_ind_count+length(vessel_graph.link.cc_ind{u})-1;
end
% edges_ind_count = edges_ind_count+length(vessel_graph.link.pos_ind);
% idx = find(link_cc_ind == 0);
% for u = 1:length(idx)
%     
% end
edges_ind = zeros(edges_ind_count,2);
%% Assign nodes and edges. 
% Convert from UCSD graph structure to Boas Lab graph structure.
% The goal is to retrieve a standard format for nodes and edges. This
% allows us to use standard toolboxes for graph theory.
% This section has been validated.
% TODO: preallocate variable "tttt"

node_idx = 1;
edge_idx = 1;
tttt = [];
link_cc_ind = zeros(size(vessel_graph.link.cc_ind));
for u = 1:length(vessel_graph.node.cc_ind)
    nodes_ind(node_idx) = vessel_graph.node.cc_ind{u}(1);
    tttt = [tttt, node_idx];
    temp_node = nodes_ind(node_idx);
    for v = 1:length(vessel_graph.node.connected_link_label{u})
        connected_link = vessel_graph.node.connected_link_label{u}(v);
        connected_link_endnodes = [vessel_graph.link.cc_ind{connected_link}(1) vessel_graph.link.cc_ind{connected_link}(end)];
        [n1,n2,n3] = ind2sub(angio_size,temp_node);
        for w = 1:2
            [l1,l2,l3] = ind2sub(angio_size,connected_link_endnodes(w));
            d(w) = sqrt((n1-l1)^2+(n2-l2)^2+(n3-l3)^2);
        end
        [~,min_idx] = min(d);
        edges_ind(edge_idx,:) = [temp_node connected_link_endnodes(min_idx(1))];
        if link_cc_ind(connected_link) == 0
            link_length = length(vessel_graph.link.cc_ind{connected_link});
            edges_ind(edge_idx+1:edge_idx+link_length-1,1) = vessel_graph.link.cc_ind{connected_link}(1:end-1);
            edges_ind(edge_idx+1:edge_idx+link_length-1,2) = vessel_graph.link.cc_ind{connected_link}(2:end);
            edge_idx = edge_idx+link_length;
            nodes_ind(node_idx+1:node_idx+link_length) = vessel_graph.link.cc_ind{connected_link};
            isa = ismember(nodes_ind(node_idx+1:node_idx+link_length),0);
            tttt = [tttt, node_idx+1:node_idx+link_length];
            node_idx = node_idx+link_length;
            link_cc_ind(connected_link) = 1;
        else
            edge_idx = edge_idx+1;
        end
    end
    node_idx = node_idx+1;
end

idx = find(link_cc_ind == 0);
for u = 1:length(idx)
    link_length = length(vessel_graph.link.cc_ind{idx(u)});
    edges_ind(edge_idx+1:edge_idx+link_length-1,1) = vessel_graph.link.cc_ind{idx(u)}(1:end-1);
    edges_ind(edge_idx+1:edge_idx+link_length-1,2) = vessel_graph.link.cc_ind{idx(u)}(2:end);
    edge_idx = edge_idx+link_length;
    nodes_ind(node_idx+1:node_idx+link_length) = vessel_graph.link.cc_ind{idx(u)};
    isa = ismember(nodes_ind(node_idx+1:node_idx+link_length),0);
    tttt = [tttt node_idx+1:node_idx+link_length];
    node_idx = node_idx+link_length;
end

%%% Create final nodes and edges variables
[n1, n2, n3] = ind2sub(angio_size,nodes_ind);
nodes =[n1', n2', n3'];
edges = zeros(size(edges_ind));

% Find indices for edges (corresponding node index).
for u = 1:size(edges_ind,1)
    edges(u,1) = find(nodes_ind == edges_ind(u,1));
    edges(u,2) = find(nodes_ind == edges_ind(u,2));
end
Graph.nodes = nodes;
Graph.edges = edges;

%% Remove redundant edges (edge = node i connected to node i)
sameEdgeIdx = [];
for u = 1:size(Graph.edges,1)
    if Graph.edges(u,1) == Graph.edges(u,2)
        sameEdgeIdx = [sameEdgeIdx; u];
    end
end
Graph.edges(sameEdgeIdx,:) = [];

%% Save the graph
% Swap x and y to convert into standard format
temp = Graph.nodes(:,2);
Graph.nodes(:,2) = Graph.nodes(:,1);
Graph.nodes(:,1) = temp;

%%% Save graph
% Remove .mat or .tif extension
filename = filename(1:end-4);
% Append correct extension
filename = strcat(filename, '_frangi_seg.mat');
fout = strcat(dpath, filename);
save(fout,'Graph');