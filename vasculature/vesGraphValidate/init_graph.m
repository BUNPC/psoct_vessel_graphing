function [Data] = init_graph(graph_nodes_segs)
%INIT_GRAPH initializes the graph data that is normally performed in GUI.
% Run the following steps:
%   - Load the graph data
%   - Run "Verification > get segment info > Update"
%   - Run "Update branch info" (Data.Graph)
%   Note: there are bugs for the following functions (regraph, prune,
%   straighten), so for now they will be excluded from this funciton.
%   - Run "Regraph Nodes" to down sample (Data.Graph_ds)
%   - Run prune_loops and prune_segment (Data.Graph_pr)
%   - Run straighten (Data.Graph_ds_pr_st)
% INPUTS:
%      graph_nodes_segs (struct): matlab struct containing the nodes and
%                                   segments.
% OUTPUTS:
%      Data (struct): graph with metadata

%% Load the graph data
Data.Graph = graph_nodes_segs;

%% Run "Verification > get segment info > Update"
%%% remove single floating nodes before calling nodeGrps_vesSegment
node_pos = Data.Graph.nodes;
edge_node_idx = Data.Graph.edges;
% Calculate number of bifurcations at each node
nbifur = zeros(size(node_pos,1),1);
n_nodes = size(node_pos,1);
% TODO: this can be replaced with "ismember" call
% ~ismember([1:n_nodes], edge_node_idx)
for ii=1:n_nodes
   nbifur(ii) = length(find(edge_node_idx(:,1)==ii | edge_node_idx(:,2)==ii));
end
% If an element of nbifur==0, then the node was not connected to any edge.
% This indicates it was not connected to any segments.
rm_list = find(nbifur == 0);
[node_pos_reindex, edge_node_reindex] =...
    remove_reindex_nodes(rm_list, node_pos, edge_node_idx);
% Update Data with reindex nodes and edges
Data.Graph.nodes = node_pos_reindex;
Data.Graph.edges = edge_node_reindex;
Data.Graph.verifiedNodes = zeros(size(Data.Graph.nodes,1),1);
%%% Call function to compute graph properties
Data.Graph.segInfo =...
    nodeGrps_vesSegment(Data.Graph.nodes, Data.Graph.edges, Data.Graph.vox);
segments = 1:size(Data.Graph.segInfo.segEndNodes,1);
segCGrps = zeros(1,size(Data.Graph.segInfo.segEndNodes,1));% group number for all segments
grpN = 1;
while ~isempty(segments)
    cSeg = segments(1);
    completedSegments = cSeg;
    segEndNodes = Data.Graph.segInfo.segEndNodes(cSeg,:);
    cdeletednodes = [];
    while ~isempty(segEndNodes)
        cEndnode = segEndNodes(end);
        seg1 = find((Data.Graph.segInfo.segEndNodes(:,1) == cEndnode) | (Data.Graph.segInfo.segEndNodes(:,2) == cEndnode));
        cSeg = unique([cSeg ; seg1]);
        tempseg = unique(setdiff(seg1,completedSegments)); 
        cdeletednodes = [cdeletednodes segEndNodes(end)];
        segEndNodes(end) = [];
        for u = 1:length(tempseg)
            tempEndNodes = Data.Graph.segInfo.segEndNodes(tempseg(u),:);
            segEndNodes = setdiff(unique([segEndNodes; tempEndNodes']),cdeletednodes);
        end
        temparray = ismember(segments,cSeg);
        completedSegments = unique([completedSegments; cSeg]);
        idx = find(temparray == 1);
        segCGrps(cSeg) = grpN;
        segments(idx) = [];
    end
    grpN = grpN+1;
end
Data.Graph.segInfo.segCGrps = segCGrps;

%% Run "Update branch info" (graph)
nSeg = size(Data.Graph.segInfo.segEndNodes,1);
idx3 = [];            A_idx3 = [];
idx5 = [];            A_idx5 = [];
idx10 = [];           A_idx10 = [];
idx_end = cell(4,1);  A_idx_end = cell(4,1);
iall = 1; i3 = 1; i5 = 1; i10 = 1;
A_iall = 1; A_i3 = 1; A_i5 = 1; A_i10 = 1;
for u = 1:nSeg
    temp_idx = find(Data.Graph.segInfo.nodeSegN == u);
    endnodes = Data.Graph.segInfo.segEndNodes(u,:);
    %        nc1 = length(find(Data.Graph.edges(:,1) == endnodes(1) |Data.Graph.edges(:,2) == endnodes(1)));
    %        nc2 = length(find(Data.Graph.edges(:,1) == endnodes(2) |Data.Graph.edges(:,2) == endnodes(2)));
    nc1 = length(find((Data.Graph.segInfo.segEndNodes(:,1) == endnodes(1)) | (Data.Graph.segInfo.segEndNodes(:,2) == endnodes(1))));
    nc2 = length(find((Data.Graph.segInfo.segEndNodes(:,1) == endnodes(2)) | (Data.Graph.segInfo.segEndNodes(:,2) == endnodes(2))));
    if nc1 < 2 || nc2 < 2
        endnode = 1;
    else
        endnode = 0;
    end
    if 1 
        if endnode == 1
            idx_end{1}(iall) = u;
            iall = iall+1;
        end
        if length(temp_idx) <= 3
            idx3 = [idx3; u];
            if endnode == 1
                idx_end{2}(i3) = u;
                i3 = i3+1;
            end
        elseif length(temp_idx) <= 5
            idx5 = [idx5; u];
            if endnode == 1
                idx_end{3}(i5) = u;
                i5 = i5+1;
            end
        elseif length(temp_idx) <= 10
            idx10 = [idx10; u];
            if endnode == 1
                idx_end{4}(i10) = u;
                i10 = i10+1;
            end
            
        end
    end
    if endnode == 1
        A_idx_end{1}(A_iall) = u;
        A_iall = A_iall+1;
    end
    if length(temp_idx) <= 3
        A_idx3 = [A_idx3; u];
        if endnode == 1
            A_idx_end{2}(A_i3) = u;
            A_i3 = A_i3+1;
        end
    elseif length(temp_idx) <= 5
        A_idx5 = [A_idx5; u];
        if endnode == 1
            A_idx_end{3}(A_i5) = u;
            A_i5 = A_i5+1;
        end
    elseif length(temp_idx) <= 10
        A_idx10 = [A_idx10; u];
        if endnode == 1
            A_idx_end{4}(A_i10) = u;
            A_i10 = A_i10+1;
        end
    end
      
end

% Number of bifurcations at end node
endnodes1 = unique(Data.Graph.segInfo.segEndNodes(:,1));
endnodes2 = unique(Data.Graph.segInfo.segEndNodes(:,2));
unique_endnodes = unique([endnodes1;endnodes2]);
Data.Graph.endNodes = zeros(size(unique_endnodes));
nB = zeros(size(unique_endnodes));
for ii=1:length(unique_endnodes)
    nB(ii) = length(find(Data.Graph.segInfo.segEndNodes(:)==unique_endnodes(ii))); 
    Data.Graph.endNodes(ii) = unique_endnodes(ii);
end
Data.Graph.nB = nB;
Data.Graph.verifiedSegments = zeros(size(Data.Graph.segInfo.segLen));
idx = find(Data.Graph.verifiedSegments == 0);
Data.Graph.segInfo.unverifiedIdx = idx;
Data.Graph.segInfo.unverifiedIdx3 = idx3;
Data.Graph.segInfo.unverifiedIdx5 = idx5;
Data.Graph.segInfo.unverifiedIdx10 = idx10;
Data.Graph.segInfo.idx_end = idx_end;
Data.Graph.segInfo.A_Idx3 = A_idx3;
Data.Graph.segInfo.A_Idx5 = A_idx5;
Data.Graph.segInfo.A_Idx10 = A_idx10;
Data.Graph.segInfo.A_idx_end = A_idx_end;
if isfield(Data.Graph.segInfo,'segCGrps')
    grpLengths = zeros(1,length(Data.Graph.segInfo.segCGrps(:)));
    for u = 1:length(Data.Graph.segInfo.segCGrps)
        idx = find(Data.Graph.segInfo.segCGrps == u);
        grpLengths(u) = length(idx);
    end
    [~,id] = sort(grpLengths,'descend');
    segmentsGrpOrder = zeros(length(Data.Graph.segInfo.segCGrps),2);
    segmentsUnverifiedGrpOrder = zeros(length(find(Data.Graph.verifiedSegments == 0)),2);
    sidx = 1;
    usidx = 1;
    for u = 1:length(id)
        idx = find(Data.Graph.segInfo.segCGrps == id(u));
        tempidx = find( Data.Graph.verifiedSegments == 0 );
        uidx = intersect(idx,tempidx);
        segmentsGrpOrder(sidx:sidx+length(idx)-1,1) = idx;
        segmentsGrpOrder(sidx:sidx+length(idx)-1,2) = u;
        segmentsUnverifiedGrpOrder(usidx:usidx+length(uidx)-1,1) = uidx;
        segmentsUnverifiedGrpOrder(usidx:usidx+length(uidx)-1,2) = u;
        sidx = sidx+length(idx);
        usidx = usidx+length(uidx);
    end
    % After trashing deleted segments and then running "Update Branch
    % Info," this next line sets the last row of the array
    % segmentsGrpOrder to zeros.
    Data.Graph.segInfo.segmentsGrpOrder = segmentsGrpOrder;
    Data.Graph.segInfo.segmentsUnverifiedGrpOrder = segmentsUnverifiedGrpOrder;
end


end