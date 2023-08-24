function flag = nodesConnected( nodeEdges, newNode, lastNode, nDis )

nDis2 = 0;
lstEdges = find(nodeEdges(:,1)==newNode | nodeEdges(:,2)==newNode);
lstDis = nDis2 * ones(length(lstEdges),1);
lstNode = newNode *  ones(length(lstEdges),1);

flag = 0;
flag2 = 0;
%flag3 = 0;

n0 = newNode;
while ~isempty(lstEdges)% | flag3
    if nDis2<nDis & ~flag2
%        flag3 = 0;
        
        nDis2 = nDis2 + 1;
        n1 = setdiff(nodeEdges(lstEdges(end),:),n0); 
        if n1==lastNode
            flag = 1;
            return
        end
        if nDis2<nDis
            lst = setdiff(find(nodeEdges(:,1)==n1 | nodeEdges(:,2)==n1),lstEdges(end));
            if isempty(lst)
                flag2 = 1;
            end
        else
            lst = [];
        end
        lstEdges = [lstEdges(1:end-1); lst];
        lstDis = [lstDis(1:end-1); nDis2*ones(length(lst),1)];
        lstNode = [lstNode(1:end-1); n1 *  ones(length(lst),1)];
        n0 = n1;
    else
        flag2 = 0;

        if ~isempty(lstNode)
            n0 = lstNode(end);
%            lstNode = lstNode(1:end-1);
            nDis2 = lstDis(end);
%            lstDis = lstDis(1:end-1);
%            n1 = setdiff(nodeEdges(lstEdges(end),:),n0);
%            lstEdges = lstEdges(1:end-1);
%            n0 = n1;
%            flag3 = 1;
        end
    end
end

        