function plotTreeSeg_select( edgeSel )
global im

E = im.nodeEdges;
V = im.nodePos;

hlSel = [];
for ii=1:length(im.Tree) 
    jj = im.Tree(ii);
    hl=plot3( V(E(jj,:),2), V(E(jj,:),1), V(E(jj,:),3), '.-' );
    set(hl,'ButtonDownFcn',sprintf('plotTreeSeg_select(%d)',jj));
    if jj == edgeSel
        set(hl,'color','r');
        hlSel = hl;
    end
end

ch = menu('Action?','Delete edge','Recenter','Collapse','Cancel');

if ch==1
    % DELETE EDGE
    edgeFlag = ones(size(E,1),1);
    edgeFlag(edgeSel) = 0;
    lstEdges = find(edgeFlag==1);
    im.nodeEdges = E(lstEdges,:);
    im.edgeFlag = im.edgeFlag(lstEdges);

    % check for and remove abandoned nodes
    nodeDel = 0;
    if length(find(E==E(edgeSel,1)))==1
        nodeDel = E(edgeSel,1);
    elseif length(find(E==E(edgeSel,2)))==1
        nodeDel = E(edgeSel,2);
    end
    
    if nodeDel ~= 0
        edgeFlag = ones(size(V,1),1);
        edgeFlag(nodeDel) = 0;
        
        nNodes = 0;
        nodePosTmp = [];
        nodeMap = zeros(size(im.nodePos,1),1);
        for iN = 1:size(im.nodePos,1)
            if edgeFlag(iN)==1
                nNodes = nNodes + 1;
                nodePosTmp(nNodes,:) = im.nodePos(iN,:);
                nodeDiamTmp(nNodes) = im.nodeDiam(iN);
                nodeMap(iN) = nNodes;
            end
        end
        im.nodePos = nodePosTmp;
        im.nodeDiam = nodeDiamTmp;
        im.nodeEdges = nodeMap(im.nodeEdges);
        im.nBflag = 1;
    end
    
    % remove redundant edges
    % point edges
    nodeEdges = im.nodeEdges;
    nodeEdges = nodeEdges(find(nodeEdges(:,1)~=nodeEdges(:,2)),:);
    % redundant edges
    sE = cell(size(nodeEdges,1),1);
    for ii=1:length(nodeEdges)
        if nodeEdges(ii,1)<nodeEdges(ii,2)
            sE{ii} = sprintf('%05d%05d',nodeEdges(ii,1),nodeEdges(ii,2));
        else
            sE{ii} = sprintf('%05d%05d',nodeEdges(ii,2),nodeEdges(ii,1));
        end
    end
    [b,i,j]=unique(sE);
    im.nodeEdges = nodeEdges(sort(i),:);
    
    im.edgeFlag = zeros(size(im.nodeEdges,1),1);

    

    I = im.III;
    lst = find(I>=65);
    I(lst) = 0;

    [pMin,pMax,im.Tree] = plotTreeSeg(im.nodeEdges,im.nodePos,im.nodeSelected,im.TreeNedges,[1 im.TreeThresh],I);
    
    hVesGraph = findobj('tag','vesselGraph');
    if ~isempty(hVesGraph)
        h = guidata(hVesGraph);
        set(h.pushbuttonImageGraph,'enable','on');
    end

elseif ch==2
    % RECENTER TREE    
    I = im.III;
    lst = find(I>=65);
    I(lst) = 0;

    im.nodeSelected = E(edgeSel,1);
    [pMin,pMax,im.Tree] = plotTreeSeg(im.nodeEdges,im.nodePos,im.nodeSelected,im.TreeNedges,[1 im.TreeThresh],I);

elseif ch==3
    % COLLAPSE NODES
    % collapse around node with most branches, if equal then choose 1
    pts = E(edgeSel,:);
    if length(find(E==pts(2))) > length(find(E==pts(1)))
        pt = pts(2);
    else
        pt = pts(1);
    end
    
    nNodes = size(V,1);
    rho = sum( (V - ones(nNodes,1)*V(pt,:) ).^2,2 ).^0.5; 
    lst = find(rho<im.nodeDiam(pt));
    
    for ii=1:length(lst)
        figure(10)
        plot3(V(lst(ii),2),V(lst(ii),1),V(lst(ii),3),'r*')
    end
    plot3(V(pt,2),V(pt,1),V(pt,3),'g*');
    
    ch = menu('Collapse these points?','Yes','No');
    if ch==1
        % collapse edges to selected node
        for ii=1:length(lst)
            if lst(ii)~=pt
                lst2 = find(E==lst(ii));
                E(lst2) = pt;
            end
        end
        % remove abandoned nodes
        edgeFlag = ones(size(V,1),1);
        edgeFlag(lst) = 0;
        edgeFlag(pt) = 1;
        nNodes = 0;
        nodePosTmp = [];
        nodeMap = zeros(size(V,1),1);
        for iN = 1:size(V,1)
            if edgeFlag(iN)==1
                nNodes = nNodes + 1;
                nodePosTmp(nNodes,:) = V(iN,:);
                nodeDiamTmp(nNodes) = im.nodeDiam(iN);
                nodeMap(iN) = nNodes;
                if iN==pt
                    im.nodeSelected = nNodes;
                end
            end
        end
        nodeEdges = nodeMap(E);
        im.nodePos = nodePosTmp;
        im.nodeDiam = nodeDiamTmp;

        % remove redundant edges
        % point edges
        nodeEdges = nodeEdges(find(nodeEdges(:,1)~=nodeEdges(:,2)),:);
        % redundant edges
        sE = cell(size(nodeEdges,1),1);
        for ii=1:length(nodeEdges)
            if nodeEdges(ii,1)<nodeEdges(ii,2)
                sE{ii} = sprintf('%05d%05d',nodeEdges(ii,1),nodeEdges(ii,2));
            else
                sE{ii} = sprintf('%05d%05d',nodeEdges(ii,2),nodeEdges(ii,1));
            end
        end
        [b,i,j]=unique(sE);
        im.nodeEdges = nodeEdges(sort(i),:);        

        im.edgeFlag = zeros(size(im.nodeEdges,1),1);

    end
    
    % RECENTER TREE    
    I = im.III;
    lst = find(I>=65);
    I(lst) = 0;

    [pMin,pMax,im.Tree] = plotTreeSeg(im.nodeEdges,im.nodePos,im.nodeSelected,im.TreeNedges,[1 im.TreeThresh],I);
    
    hVesGraph = findobj('tag','vesselGraph');
    if ~isempty(hVesGraph)
        h = guidata(hVesGraph);
        set(h.pushbuttonImageGraph,'enable','on');
    end

else
    if ~isempty(hlSel)
        set(hlSel,'color','b');
    end
end
