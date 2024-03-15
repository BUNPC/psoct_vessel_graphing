function [sparsity, cnodes] = graph_sparsity(edges)
%GRAPH_SPARSITY determine sparsity of cycles (loops) in graph
% The purpose is to determine if the cycles have been fully down sampled.
% In this case, # nodes in the cycle == # edges in cycle. This is to
% determine whether the cycle is fully downsampled.
%
% OUTLINE:
%   - generate graph from edges
%   - identify all cycles in entire graph (allcycles() command)
%   - iterate over each cycle (cnodes is a cell array of nodes in loop)
%       - calculate degree for each node in loop
%       - if all nodes are greater than second degree
%           cycle is sparse
%       - else
%           cycle is not sparse
%
% INPUTS:
%   nodes (Nx2 array): node positions
%   edges (Nx3 array): array of edges [n1, n2]
% OUTPUTS:
%   sparsity (bool matrix):
%           1 for sparse loops
%           0 for dense loops

%% Initalization
% Convert edges to graph
g = graph(edges(:,1), edges(:,2));

% Generate arrays of node and edge indices (cnodes, cedges) for each cycle
[cnodes, ~] = allcycles(g,"MaxNumCycles",500);

% Generate bool matrix for sparsity
sparsity = zeros(length(cnodes),1);

%% Identify sparsity of each cycle

% Iterate over cell array of cycles
for ii=1:length(cnodes)
    % extract # nodes
    n = cnodes{ii,:};
    
    % Find degree of each node in cycle
    d = degree(g, n);

    %%% If degree of all nodes > 2, then cycle is sparse
    if all(d>2)
        sparsity(ii) = true;
    else
        sparsity(ii) = false;
    end
end

end