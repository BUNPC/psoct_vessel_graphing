function seg_graph_init(seg, vox_dim, fullpath, fname_seg)
% Initialize graph from segmentation
% INPUTS:
%   seg (mat): segmentation matrix
%   vox_dim (array): 3-element array of voxel dimensions (microns) [x,y,z]
%   fullpath (string): output directory for graph
%   fname_seg (string): filename of segmentation

%% Convert segmentation to graph (nodes and edges)
graph = seg_to_graph(seg, vox_dim);

%% Remove loops from graph
% Move to mean minimum voxel intensity
v_min = 0.99;

% Initail search radius in down sample function (voxels)
delta = 6;

% # iterations for mv2mean function in for-loop iteration in rm_loops
mv_iter = 1;

% Call function to remove loops
[nodes_rm, edges_rm] = rm_loops(graph.nodes, graph.edges, seg, delta, v_min, mv_iter);

% Update graph with nodes and edges
graph.nodes = nodes_rm;
graph.edges = edges_rm;

%% Initialize graph metadata (Graph.Data)
[Data] = init_graph(graph);

%%% Append "angio" data (segmentation matrix)
% Rearrange [x,y,z] --> [z,x,y]. This is the standard in
% the graph validation GUI.
angio = permute(seg, [3,1,2]);
Data.angio = angio;

% Create new filename for graph and add .MAT extension
fname_graph = strcat(fname_seg, '_graph_data.mat');
fout = fullfile(fullpath, fname_graph);
save(fout,'Data', '-v7.3');
end