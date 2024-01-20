function im = mv_to_mean(im, v_min, varargin)
%mv_to_mean Move nodes towards mean of neighboring nodes (collapse loops)
% Outline:
%   - iterate over every node index in the graph structure.
%   - "current node index" is the ii'th node index in the list of total node indices.
%   - posO = location of current node
%   - posC = mean location of connected nodes
%   - posN = new position as average of norm(posC - posO)
%   - verify that posN is within v_min (voxel intensity)
%
%   INPUTS:
%       im.nodes ([n,3] array): node locations
%       im.edges ([m,2] array): edges connecting each node
%       im.angio (double matrix): PS-OCT intensity volume (vessels are
%               bright)
%       v_min (double): minimum voxel intensity threshold. The new voxel
%               position will only be reassigned if the voxel intensity of
%               the new node position is >= v_min.
%       varargin: node indices belonging to cycle, specified as [1, N]
%   OUTPUTS:
%       n ([n,3] array): node locations
%       e ([m,2] array): edges connecting each node


%% Assign local variables
[x,y,z] = size(im.angio);
nodes = im.nodes;

%%% Check if using subset of nodes or all nodes
if ~isempty(varargin)
    nLst = varargin{1};
else
    nLst = 1:size(nodes,1);
end
hwait = waitbar(0,'Moving nodes towards mean of neighboring nodes');
dbstop if error
%% Iterate over nodes
for ii = 1:length(nLst)
    % Retrieve node index
    nidx = nLst(ii);
    waitbar(ii/length(nLst),hwait)
    
    %%%% Find edges containing current node index
    eLst = find(im.edges(:,1)==nidx | im.edges(:,2)==nidx);
    
    %%% Find indices of nodes connected to current node index
    % For the edges connected to the current node index, find the indices
    % of other nodes connected to the same edge "unique(im.edges(eLst,:))"
    % Use setdiff to exclude the index of the current node index (nidx).
    nLst2 = setdiff(unique(im.edges(eLst,:)), nidx);

    %%% Iterate over indices of connected nodes
    if length(nLst2)>1
        % Set null position to currrent node index or [1,1,1]
        pos0 = max(im.nodes(nidx,:),1);
        % Set C position to mean of connected nodes
        posC = mean(im.nodes(nLst2,:),1);
        % New position = (null position + difference) / Euclidean distance
        posN = pos0 + (posC-pos0) / max(norm(posC-pos0),1);
        % If the graph node is outside of segmentation, then pass
        if round(posN(1)) > x || round(posN(2)) > y || round(posN(3)) > z
            return
        % If new position is within vessel (>=vox. intensity threshold)
        % Then reassign node position to new position.
        elseif im.angio(round(posN(1)),round(posN(2)),round(posN(3)))>=v_min
            nodes(nidx,:) = posN;
        end
    end
end
im.nodes = nodes;
close(hwait)
end