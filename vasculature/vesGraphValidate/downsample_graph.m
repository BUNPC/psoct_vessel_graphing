function [nodes, edges] =...
    downsample_graph(nodes, edges, vox_xy, vox_z)
%%% Down sample 
% INPUTS:
%       nodes (matrix): Nodes prior to downsample
%       edges (matrix): edges prior to downsample
%       validatedNodes (array): validated node indices
%       vox_xy (double): voxel dimensions (x,y)
%       vox_z (double):  voxel dimension (z)
% OUTPUTS:
%       nodes (matrix): 
%       edges (matrix): 
%       validatedNodes (array): validated node indices
%       validatedEdges (array): validated edges
%{
TODO:
- rewrite this function so that the outer for-loop iterates over the
segments. Within a segment, iterate over the nodes and compare the location
of nodes within the boundaries [hxy, hxy, hz].
%}



%% Initialization



fprintf('Regraph reduced %d nodes to %d\n',nNodes,size(nodes,1))

