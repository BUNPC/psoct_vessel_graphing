function seg_to_graph(dpath, filename)
%seg_to_graph Convert the vessel segments to a graph (nodes + vertices)
% This function is the second step in the vessel segmentation pipeline.
% This function can take either a .MAT or .TIF
%
%%% INPUTS:
%       dpath (str): absolute path to the data folder
%       fname (str): name of the segmentation file (.MAT or .TIFF)
%%% OUTPUTS:
%       This function does not return anything. Instead, it saves a struct
%       to the same directory that contains the segmentation files. The
%       name of the struct will be "fname_frangi_seg"


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
newdir = mydir(1:idcs(end-1));
addpath(genpath(newdir));

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

%%% This function output is unused. It is leftover from original code.
% Compute Euclidean distance transform of inverse of binary vessel mask
% vessel_mask_dt = bwdist(~vessel_mask);

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

edges_ind = zeros(edges_ind_count,2);
%% Assign nodes and edges. 
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

%% TODO: determine purpose of this section
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

%% TODO: determine purpose of this section

[n1, n2, n3] = ind2sub(angio_size,nodes_ind);
nodes = [n1', n2', n3'];
edges = zeros(size(edges_ind));

%% TODO: find faster search method
for u = 1:size(edges_ind,1)
    edges(u,1) = find(nodes_ind == edges_ind(u,1));
    edges(u,2) = find(nodes_ind == edges_ind(u,2));
end

%% Remove redundant edges from graph
Graph.nodes = nodes;
Graph.edges = edges;

% Create array for storing indices of edges with the same edge.
nedge = size(Graph.edges,1);
same_edge_idx = zeros(nedge,1);
j = 1;      % index for position in same_edge_idx array

% Iterate over all edges
for ii = 1:nedge
    if Graph.edges(ii,1) == Graph.edges(ii,2)
        same_edge_idx(j) = ii;
        j = j + 1;
    end
end

% Remove zero entries in same_edge_idx. Leave just nonzero indices.
same_edge_idx(same_edge_idx==0) = [];

% Delete the edges with a shared edge
Graph.edges(same_edge_idx,:) = [];

%% Save the graph

% Update graph
temp = Graph.nodes(:,2);
Graph.nodes(:,2) = Graph.nodes(:,1);
Graph.nodes(:,1) = temp;

%%% Save graph
% Remove .mat or .tif extension
filename = filename(1:end-4);
% Create new filename for graph and add .MAT extension
filename = strcat(filename, '_frangi_seg.mat');
fout = strcat(dpath, filename);
save(fout,'Graph');

end