function im = mv_to_mean(im, Ithresh)
%% Move nodes towards mean of neighboring nodes
% This script will iterate over every node index in the graph structure.
% The comments reference "current node index," which is referring to
% the ii'th node index in the list of total node indices.

% Assign local variables
nodes = im.nodes;
nLst = 1:size(nodes,1);
hwait = waitbar(0,'Moving nodes towards mean of neighboring nodes');

%% Iterate over nodes
for ii = 1:size(nodes,1)
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
        % If new position is within vessel (>=vox. intensity threshold)
        % Then reassign node position to new position.
        if im.angio(round(posN(1)),round(posN(2)),round(posN(3)))>=Ithresh
            nodes(nidx,:) = posN;
        end
    end
end
im.nodes = nodes;
close(hwait)
end