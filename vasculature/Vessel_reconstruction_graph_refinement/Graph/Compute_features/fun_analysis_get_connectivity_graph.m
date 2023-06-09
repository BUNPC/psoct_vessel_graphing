function graph_str = fun_analysis_get_connectivity_graph(vessel_graph)
% fun_analysis_get_connectivity_graph construct MATLAB graph structures and
% adjacent matrices of weighted, unweighted and link label from the input
% graph generated by fun_skeleton_to_graph
% Only the links without endpoint(s) are used for the graph construction. 
% Input: 
%   vessel_graph: structure generated by fun_skeleton_to_grah
% Output: 
%   graph_str: structure with fields: 
%       graph_uw: unweighted graph structure
%       grpah_w: weighted graph structure, weights of link is defined as
%       the physicsal length of the link segments
%       A_weighted: weighted adjacent matrix derived from graph_w
%       A_link_label: symmetric sparse matrix whose (i,j) element is the link label
%       between node i and node j. Multiedges has be removed. 
%       used_link_label: label of the links in input_graph that are used
%       for graph construction. 
%       bilink_loop: structure generated by fun_graph_get_bilink_loops,
%       with additional link length data. 
% Implemented by Xiang Ji on 02/24/2019
% To do: select one of the link in the bilink-loop by some metric. Now it
% is selected by label order, which is not very reasonable. 
% Note: this function only work for undirectly graph!
% Only construct graph using links without endpoints

% Modified by Xiang Ji on 03/27/2019
% 1. Compute the length of the link as edge weight inside the function -
% remove dependence on the vessel_graph.link.features.length
% 
no_endpoint_Q = all(vessel_graph.link.connected_node_label, 2);
% Compute the length of the link here
tmp_sub = cellfun(@(x) fun_ind2sub(vessel_graph.num.mask_size, x), vessel_graph.link.cc_ind, 'UniformOutput', false);
all_edge_weight =  cellfun(@fun_graph_sub_to_length, tmp_sub);

no_endpoint_label = find(no_endpoint_Q);
% Assume vascular network does not have self-loop - should be reasonable
self_loop_str = fun_graph_get_self_loops(vessel_graph);
if ~isempty(self_loop_str.link_label)
%     assert(isempty(self_loop_str.link_label), 'Exist self-loop');
%     warning('Exist self-loop. Removed it before constructing the grpah');
    no_endpoint_label = setdiff(no_endpoint_label, self_loop_str.link_label, 'stable');
end
graph_str.used_link_label = no_endpoint_label;

valid_connected_node_label = vessel_graph.link.connected_node_label(no_endpoint_label, :);
graph_str.used_node_label_pair = valid_connected_node_label;

valid_edge_weight = all_edge_weight(no_endpoint_label);
graph_str.used_link_length = valid_edge_weight;

% Constructed undirected graph 
graph_str.graph_uw = graph(valid_connected_node_label(:,1), valid_connected_node_label(:,2));
graph_str.graph_w = graph(valid_connected_node_label(:,1), valid_connected_node_label(:,2), valid_edge_weight);
% Construct the adjacent matrix for links that are not in the self-loop nor
% in the bi-link loop 

[tmp_node_label, tmp_idx, ~] = unique(sort(valid_connected_node_label, 2, 'ascend'), 'rows', 'stable');

tmp_link_label = no_endpoint_label(tmp_idx);
tmp_sub_1 = tmp_node_label(:);
tmp_sub_2 = [tmp_node_label(:,2); tmp_node_label(:,1)];
tmp_val = [tmp_link_label; tmp_link_label];
tmp_weight = [valid_edge_weight(tmp_idx);valid_edge_weight(tmp_idx)];

graph_str.A_link_label = sparse(tmp_sub_1, tmp_sub_2, tmp_val, vessel_graph.node.num_cc, ...
    vessel_graph.node.num_cc);

graph_str.A_weighted = sparse(tmp_sub_1, tmp_sub_2, tmp_weight, vessel_graph.node.num_cc, ...
    vessel_graph.node.num_cc);

graph_str.Adj_mat_size = [vessel_graph.node.num_cc, vessel_graph.node.num_cc];

% Deal with loops consists of two link - i.e. multiedge case 
bilink_loop_str = fun_graph_get_bilink_loops(vessel_graph);
% Add length information to the bilink loop structure:
bilink_loop_str.link_length_pair = all_edge_weight(bilink_loop_str.link_label_pair);
tmp_num_iter = numel(bilink_loop_str.link_label_pair_ge_3);
bilink_loop_str.link_length_pair_ge_3 = cell(tmp_num_iter, 1);
for iter = 1 : tmp_num_iter
    bilink_loop_str.link_length_pair_ge_3{iter} = all_edge_weight(...
        bilink_loop_str.link_label_pair_ge_3{iter});
end
graph_str.bilink_loop = bilink_loop_str;
end