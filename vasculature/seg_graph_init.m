function seg_graph_init(seg, vox_dim, fullpath, fname_seg, viz, rmloop_bool)
% Initialize graph from segmentation
% INPUTS:
%       seg (mat): segmentation matrix
%       vox_dim (array): 3-element array of voxel dimensions (microns) [x,y,z]
%       fullpath (string): output directory for graph
%       fname_seg (string): filename of segmentation
%       viz (bool): true = display debugging graph figures.
%       rmloop_bool (bool): true = remove loops from graph

%% Convert segmentation to graph (nodes and edges)
graph = seg_to_graph(seg, vox_dim);

%% Remove loops from graph
% Move to mean minimum voxel intensity
v_min = 0.99;

% Initail search radius in down sample function (voxels)
delta = 6;

% # iterations for mv2mean function in for-loop iteration in rm_loops
mv_iter = 1;

if rmloop_bool
    % Call function to remove loops
    [nodes_rm, edges_rm] = rm_loops(graph.nodes, graph.edges, seg, delta, ...
                                    v_min, mv_iter, viz);
    % Update graph with nodes and edges
    graph.nodes = nodes_rm;
    graph.edges = edges_rm;
end

%% Initialize graph metadata (Graph.Data)
[Data] = init_graph(graph);

%%% Append "angio" data (segmentation matrix)
% Note that the vesGraphValidate GUI expects the segmentation vector in a
% permuted format: [x,y,z] --> [z,x,y]. For continuity, this script does
% not permute the segmentation. If this is needed later, then the code is
% on the next line:
% seg = permute(seg, [3,1,2]);
Data.angio = seg;

%%% Create new filename for graph and add .MAT extension
if rmloop_bool
    fname_graph = strcat(fname_seg, 'rmloop_graph_data.mat');
else
    fname_graph = strcat(fname_seg, 'graph_data.mat');
end
fout = fullfile(fullpath, fname_graph);
save(fout,'Data', '-v7.3');
end