%% Remove looped segments
%{
This script uses a nested struct "Graph" which contains metadata about
the Graph struct. The Graph struct is created in the script SegtoGraph.m.

TODO:
- Determine meaning of A_idx_end{1,1}
- Determine logic of last two steps (new edges and new nodes)
%}


function [Graph] = prune_loops(Graph)

endnodes = unique(Graph.segInfo.segEndNodes(:));
numLoops = 0;
numLoopsL2 = 0; numLoops2T4 = 0; numLoopsG4 = 0;
nodeLoop = zeros( size(Graph.nodes,1), 1);
loops = [];

for uu = 1:length(endnodes)
    start_node = endnodes(uu);
    dist_cutoff = 150;   % parameter could be tuned
    dist_nodes = 0;
    edges_processed = [];
    nodes_tobeprocessed = start_node;
    nodes_processed = [];
    currentnode = start_node;
    while ~isempty(nodes_tobeprocessed) 
        while dist_nodes <= dist_cutoff
            nodes_processed = unique([nodes_processed; currentnode]);
            edges_connected = find(Graph.edges(:,1) == currentnode | Graph.edges(:,2) == currentnode);
            edges_processed = unique([edges_processed; edges_connected(:)]);
            connected_nodes = [];
            for u = 1:length(edges_connected)
                temp = Graph.edges(edges_connected(u),:);
                connected_nodes = [connected_nodes; temp(:)];
            end
            nodes_tobeprocessed = setdiff(unique([nodes_tobeprocessed;connected_nodes]),nodes_processed);
            if ~isempty(nodes_tobeprocessed) 
                currentnode = nodes_tobeprocessed(end);
            else
                break;
            end
            dist_nodes = sqrt(sum(Graph.nodes(start_node,:) - Graph.nodes(currentnode,:)).^2);
%             if ~(dist_nodes <= dist_cutoff)
%                 nodes_tobeprocessed = nodes_tobeprocessed(1:end-1);
%                 break;
%             end
        end
        nodes_processed = unique([nodes_processed; currentnode]);
        nodes_tobeprocessed = nodes_tobeprocessed(1:end-1);
        if ~isempty(nodes_tobeprocessed) 
            currentnode = nodes_tobeprocessed(end);
        end
    end  
    edges_new = [];
    for u = 1:length(nodes_processed)
        edges_new = [edges_new; find(Graph.edges(:,1) == nodes_processed(u) | Graph.edges(:,2) == nodes_processed(u))];
    end
    edges_new = unique(edges_new);
    
    %% Check if any edges detected. If not, skip this logic    
    if ~isempty(edges_new)
        esub = Graph.edges(edges_new,:);        
        v = unique(esub(:));
        esub2 = esub;
        c(uu).nodes = nodes_processed;
        for iv = 1:length(v)
            esub2(find(esub2==v(iv))) = iv;
        end
        disp(uu)
        c(uu).ind = grCycleBasis(esub2);
    
        countLoop = [];
        for iLoop = 1:size(c(uu).ind,2)
            countLoop(iLoop) = 0;
            loopNode = [];
            lstEdges = find(c(uu).ind(:,iLoop)>0);
            for iE = 1:length(lstEdges)
                if nodeLoop( esub(lstEdges(iE),1) ) == 0
                    nodeLoop( esub(lstEdges(iE),1) ) = numLoops + iLoop;
                    countLoop(iLoop) = 1;
                    loopNode = esub(lstEdges(iE),1);
                end
                if nodeLoop( esub(lstEdges(iE),2) ) == 0
                    nodeLoop( esub(lstEdges(iE),2) ) = numLoops + iLoop;
                    countLoop(iLoop) = 1;
                    loopNode = esub(lstEdges(iE),2);
                end
            end
            if ~isempty(loopNode)
                loops = [loops; loopNode];
                pt = Graph.nodes(loopNode,:);
                if pt(3) <= 100
                    numLoopsL2 = numLoopsL2+1;
                elseif pt(3) > 100 && pt(3) <= 200
                    numLoops2T4 = numLoops2T4+1;
                else
                    numLoopsG4 = numLoopsG4+1;
                end
            end
        end
        if size(c(uu).ind,2)>0
            numLoops = numLoops + length(find(countLoop>0));
        end

    end
end

Graph.segInfo.loops = loops;
disp([num2str(length(loops)) ' loops are detected in the graph.']);

%%  What is purpose of this?
segmentstodelete = [];
for v = 1:length(Graph.segInfo.loops)   
    nodeno = Graph.segInfo.loops(v);
    if ~isempty(nodeno)
        segmentstodelete = [ segmentstodelete; Graph.segInfo.nodeSegN(nodeno)];
    end
end
nodes = Graph.nodes;
edges = Graph.edges;

%%  What is purpose of this?
lstRemove=[];
for i=1:length(segmentstodelete)
    idx =  find(Graph.segInfo.nodeSegN == segmentstodelete(i));
    endnodes = Graph.segInfo.segEndNodes(segmentstodelete(i),:);
        endnodes = endnodes(:);
        segs1 = find(Graph.segInfo.segEndNodes(:,1) == endnodes(1) | Graph.segInfo.segEndNodes(:,2) == endnodes(1));
        segs2 = find(Graph.segInfo.segEndNodes(:,1) == endnodes(2) | Graph.segInfo.segEndNodes(:,2) == endnodes(2));
        if length(segs1) > 1
            tsegs1 = setdiff(segs1,segmentstodelete);
            if ~isempty(tsegs1)
                idx = setdiff(idx,endnodes(1));
            end
        end
        if length(segs2) > 1
            tsegs2 = setdiff(segs2,segmentstodelete);
            if ~isempty(tsegs2)
                idx = setdiff(idx,endnodes(2));
            end
        end
        lstRemove = [lstRemove; idx];
end
nNodes = size(nodes,1);

%% Remove nodes
% The variable lstRemove containes the nodes to be removed
% Some of these elements exceed the number of nodes in "map"

% Remove elements from lstRemove exceeding the length of "map"
rm_list = lstRemove(lstRemove <= length(map));

% Remove nodes
map = (1:nNodes)';
map( rm_list ) = [];
mapTemp = (1:length(map))';
nodeMap = zeros(nNodes,1);
nodeMap(map) = mapTemp;

edgesNew = nodeMap(edges);
nodesNew = nodes;
nodesNew(rm_list,:) = [];

zero_idx = find(edgesNew(:,1) == 0 | edgesNew(:,2)==0);
edgesNew(zero_idx,:) = [];

%% Save Data
Graph.nodes = nodesNew;
Graph.edges = edgesNew;
Graph.verifiedNodes = zeros(size(Graph.nodes,1),1);
disp('Loops removed from the graph.');
save('volume_nor_inverted_masked_sigma1_segpruned_looppruned.mat', 'Data', '-v7.3');

end

