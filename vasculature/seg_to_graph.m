function [graph] = seg_to_graph(angio, vox_dim)
%seg_to_graph Convert the vessel segments to a graph (nodes + vertices)
%%% INPUTS:
%       angio (matrix): segmented volume
%       vox_dim (array): voxel dimensions (x, y, z) (micron)
%%% OUTPUTS:
%       Graph (struct): nodes and edges for graph of segmentation


%% Create the local vessel_mask variable
vessel_mask = logical(angio);

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

% find number of edges to preallocate matrix
edges_ind_count = 0;
for u = 1:length(vessel_graph.node.connected_link_label)
    edges_ind_count =...
        edges_ind_count + length(vessel_graph.node.connected_link_label{u});
end
for  u = 1:length(vessel_graph.link.cc_ind)
    edges_ind_count =...
        edges_ind_count + length(vessel_graph.link.cc_ind{u})-1;
end
% Preallocate matrix for storing edge indice
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
% Initialize graph struct
graph.nodes = nodes;
graph.edges = edges;
% Add dimensions of voxel (microns)
graph.vox = vox_dim;

% Create array for storing indices of edges with the same edge.
nedge = size(graph.edges,1);
same_edge_idx = zeros(nedge,1);
j = 1;      % index for position in same_edge_idx array

% Iterate over all edges
for ii = 1:nedge
    if graph.edges(ii,1) == graph.edges(ii,2)
        same_edge_idx(j) = ii;
        j = j + 1;
    end
end

% Remove zero entries in same_edge_idx. Leave just nonzero indices.
same_edge_idx(same_edge_idx==0) = [];

% Delete the edges with a shared edge
graph.edges(same_edge_idx,:) = [];

%% Update graph
temp = graph.nodes(:,2);
graph.nodes(:,2) = graph.nodes(:,1);
graph.nodes(:,1) = temp;

end