function im = nodeGrps_vesSegment(nodePos, nodeEdges)
% Compute graph node and edge properties from a matrix of node positions in
% (x,y,z) coordinates and a matrix of the node endpoints for each segment.
% NOTE: this 
%
% INPUTS:
%   nodePos (double mat): [x,3] matrix of the (x,y,z) coordinates for
%                           each node in the graph. total number nodes = x
%   nodeEdges (double mat): [y,2] matrix. Each row corresponds to an edge
%                             in the graph, and the two values equal the
%                             node index that the edge connects. For
%                             example, if the first entry is [1,2], then
%                             this indicates the first edge is connecting
%                             node 1 and node 2.
% OUTPUTS:
%   im.nodeGrp (double vector): node group index for each node. the values
%                               range from [1, number of nodes]
%   im.nodeSegN (double vector): segment number corresponding to node
%   im.segNedges (double vector): # edges contained within each segment.
%                                 For example, a segment can connect two end
%                                 nodes and there can be multiple nodes and
%                                 edges in between.
%   im.segLen (double vector): length of vector in pixels?
%   im.segLen_um (double vector): length of vector in microns
%   im.edgeSegN (double vector): Segment number for each edge
%   im.segEndNodes (double mat): The end node indices for each segment.
%   im.segPos (double vector): Average (x,y,z) position for each segment.
%                              The average is calculated from the average
%                              of the the two end nodes (x,y,z) coordinates

%%%%%%%%%%
% TO DO:
% LOW PRIORITY:
% correct eLen (segLength) to be in um (non-cubical voxels problem) -  easy
% correct segDiam to be in um (non-cubical voxels problem) - more
% complicated - we already have calculated diameters in imView in voxels
%%%%%%%%%%


%% nB = number of bifurcations at each segment
nB = zeros(size(nodePos,1),1);
nN = size(nodePos,1);
for ii=1:nN
   nB(ii) = length(find(nodeEdges(:,1)==ii | nodeEdges(:,2)==ii));
end

%% Logic for removing single floating nodes.
% This is performed in vesGraphValidate, so it should not be necessary.
%{
lstRemove = find(nB == 0);

% Call "remove_reindex_nodes" here
nNodes = size(nodePos,1);
map = (1:nNodes)';
map(lstRemove) = [];
mapTemp = (1:length(map))';
nodeMap = zeros(nNodes,1);
nodeMap(map) = mapTemp;

edgesNew = nodeMap(nodeEdges);
[ir,~] = find(edgesNew == 0);
edgesNew(ir,:) = [];

nodePos(lstRemove,:) = [];
nodeEdges = edgesNew;

% After calling function:
% Remove disconnected nodes from list of bifurcations.
nB(lstRemove,:) = [];
% Recalculate number of nodes.
nN = size(nodePos,1);
%}

%% create nodePos_um variable, to have positions in um
% TODO:
%   - replace hard-coded hxy/hz with variable inputs to function
%   - update im.nodePos_um whenever we update im.nodePos
hxy = 2;
hz  = 2;
nodePos_um = nodePos; 
nodePos_um(:,1:2) = nodePos(:,1:2).*hxy;
nodePos_um(:,3) = nodePos(:,3).*hz;

%%% Create a list of end nodes based upon number of bifurcations at each
%%% node. Bifurcations = 2 indicates node is within a segment.
lst3p = find(nB>=3 | nB==1);

%%% Initialize matrices
edgeSegN = zeros(size(nodeEdges,1),1); % segment # for each edge
nodeSegN = zeros(size(nodePos,1),1); % segment # for each node
nodeGrp = zeros(nN,1);  % Initial array of group numbers
segNedges = [];     % number of edges contained within each segment
segEndNodes = [];   % end node indices for each segment.
nSeg = 1;   % segment index for each node
nGrp = 0;   % group index for each node

%% Iterate through list of end nodes w/ 1 or >=3 bifurcations
for ii = 1:length(lst3p)
    % Display the current iteration number to console
    sprintf('Iteration %u out of %u', ii, length(lst3p));

    % iN is the index of the index of the segment with 1 or 3 bifurcations
    iN = lst3p(ii);

    % Original Code:
    % If group number equals zero, then it has not yet been assigned.
    % Assign the variable "nGrp" to nodeGrp(iN).
    if nodeGrp(iN)==0
        nGrp = nGrp + 1;
        nodeGrp(iN) = nGrp;
    end
    
    % set nGrpN = nGrp
    nGrpN = nodeGrp(iN);

    %% The code is breaking due to logic for creating lstE
    % List of edges connected to nodes with indices matching iN. These
    % edges contain 3 bifurcations. This is one source of the zero-index
    % issue. When none of the elements in nodeEdges match these conditions,
    % the array lstE is empty. Therefore, the second for-loop is skipped
    % and nodeSegN remains 0 for this element.

    % Iterating the variable "nSeg" leads to a dimension mismatch b/w 
    % im.segEndNodes and segLen
    % squeeze(mean(reshape(nodePos(im.segEndNodes,:),[2 length(segLen) 3]),1));
    
    lstE = find(nodeEdges(:,1)==iN | nodeEdges(:,2)==iN);
    if isempty(lstE)
        nodeSegN(iN) = nSeg;
        nSeg = nSeg + 1;
    end

    % iNstart tracks the value of iN prior to entering second for loop
    % iNstart is only used in after while loop "nodeSegN(iNstart) = nSeg;"
    iNstart = iN;

    %% TODO: determine purpose of for-loop
    for jj = 1:length(lstE)
        iN = lst3p(ii);
        eIdx =lstE(jj);
        % correct this to be in um (non-cubical voxels)
        eLen = sum(diff( nodePos(squeeze(nodeEdges(eIdx,:)),:), 1, 1).^2).^0.5;
        eLen_um = sum(diff( nodePos_um(squeeze(nodeEdges(eIdx,:)),:), 1, 1).^2).^0.5;    

        % Find unique elements in nodeEdges
        iN = setdiff(unique(nodeEdges(eIdx,:)), iN);
        if nodeGrp(iN)==0
            nodeGrp(iN) = nGrpN;
        elseif nodeGrp(iN) < nGrpN
            nodeGrp(nodeGrp==nGrpN) = nodeGrp(iN);
            nGrpN = nodeGrp(iN);
        elseif nodeGrp(iN) > nGrpN
            nodeGrp(nodeGrp==nodeGrp(iN)) = nGrpN;
%            if nB(iN)<3
%                error('We should never get here');
%            end
        end
        
        % Old code that was commented.
%{
        if nodeSegN(iN)==0  % remove 6/3/09 since it appears below
            nodeSegN(iN) = nSeg; % and should resolve an issue in
                                 % nodeGrps()
            nodeSegN(iNstart) = nSeg;            
        end
%}
        
        % Original comment: re-added (8/3/09) since nodeSegN=0 indicated
        % one node from bifurcation.
        if nodeSegN(iN)==0 
            nodeSegN(iN) = nSeg;

        end

        nE = 1;
        nLst = [];

        %% TODO: determine purpose of if-statement and while-loop
        if edgeSegN(eIdx)==0
            while nB(iN)==2
                nLst(end+1) = iN;
                edgeSegN(eIdx) = nSeg;
                lstE2 = find(nodeEdges(:,1)==iN | nodeEdges(:,2)==iN);
                eIdx = setdiff(lstE2,eIdx);
                % correct this to be in um (non-cubical voxels)
                eLen = eLen + sum(diff( nodePos(squeeze(nodeEdges(eIdx,:)),:), 1, 1).^2).^0.5;
                eLen_um = eLen_um + sum(diff( nodePos_um(squeeze(nodeEdges(eIdx,:)),:), 1, 1).^2).^0.5;
                iN = setdiff(unique(nodeEdges(eIdx,:)), iN);
                if nodeGrp(iN)==0
                    nodeGrp(iN) = nGrpN;
                elseif nodeGrp(iN) < nGrpN 
                    lst = find(nodeGrp==nGrpN);
                    nodeGrp(lst) = nodeGrp(iN);
                    nGrpN = nodeGrp(iN);
                elseif nodeGrp(iN) > nGrpN
                    lst = find(nodeGrp==nodeGrp(iN));
                    nodeGrp(lst) = nGrpN;
                end
                
                if nodeSegN(iN)==0
                    nodeSegN(iN) = nSeg;
                end
                nE = nE + 1;
            end
            
            % added 6/3/09 when deleted above
            nodeSegN(iNstart) = nSeg; 
            % added 6/3/09 to resolve issue with nodeGrps
            nodeSegN(iN) = nSeg;
            edgeSegN(eIdx) = nSeg;
            segNedges(nSeg) = nE;
            segLen(nSeg) = eLen;
            segLen_um(nSeg) = eLen_um;

            kk = find(lst3p==iN);
            if ~isempty(kk) 
                segEndNodes(end+1,:) = lst3p([ii kk]);
            else
                error('We should not get here')
            end
            nSeg = nSeg + 1;
        end
    end
end
% close(hwait)

%% Remove groups with zero nodes
nGrpN = 0;
for ii=1:nGrp
    lst = find(nodeGrp==ii);
    if ~isempty(lst)
        nGrpN = nGrpN + 1;
        nodeGrp(lst) = nGrpN;
    end
end

%% Assign metadata to struct
im.nodeGrp = nodeGrp;
im.nodeSegN = nodeSegN;
im.segNedges = segNedges;
im.segLen = segLen;
im.segLen_um = segLen_um;
im.edgeSegN = edgeSegN;
im.segEndNodes = segEndNodes;
im.segPos = squeeze(mean(reshape(nodePos(im.segEndNodes,:),[2 length(segLen) 3]),1));


end