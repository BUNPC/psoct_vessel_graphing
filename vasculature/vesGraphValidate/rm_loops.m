function [h] = rm_loops(g)
%rm_loops Remove loops in graph.
%   This should be run after downsampling the graph. This function will 
%   take the nodes and edges and create a Matlab graph data structure.
%   Then, this function uses the function "allcycles" to find the loops in
%   the graph. It then removes nodes that are in the loops.
%
%   INPUTS:
%       g (graph): downsampled graph
%   OUTPUTS:
%       h (graph): loops removed



end