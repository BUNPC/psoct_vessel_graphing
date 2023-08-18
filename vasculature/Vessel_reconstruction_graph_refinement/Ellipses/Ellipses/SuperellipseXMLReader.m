function tracegraph = SuperellipseXMLReader(fname)
xDoc = xmlread(fname);
allNodes = xDoc.getElementsByTagName('Superellipse');
for k = 0:allNodes.getLength-1
    node = allNodes.item(k);
    nodeAttribs = node.getAttributes;
    n.position = zeros(1,3);
    n.ID = -1;
    n.R = zeros(3)-1;
    n.scale = -1;
    for j = 0:nodeAttribs.getLength-1
        attrib = nodeAttribs.item(j);
        tag = char(attrib.getName);
        val = char(attrib.getValue);
        if strcmp(tag, 'ID')
            n.ID = str2num(val);
        elseif strcmp(tag, 'TraceID')
            n.TraceID = str2num(val);
        elseif strcmp(tag, 'x')
            n.position(1) = str2num(val) + 1; %one is added because matlab has 1 based indexing
        elseif strcmp(tag, 'y')
            n.position(2) = str2num(val) + 1; %one is added because matlab has 1 based indexing
        elseif strcmp(tag, 'z')
            n.position(3) = str2num(val) + 1; %one is added because matlab has 1 based indexing
        elseif strcmp(tag, 'a1')            
            n.scale(1) = str2num(val);
        elseif strcmp(tag, 'a2')            
            n.scale(2) = str2num(val);
        elseif strcmp(tag, 'a3')            
            n.scale(3) = str2num(val);
        elseif strcmp(tag, 'R11')            
            n.R(1,1) = str2num(val);
        elseif strcmp(tag, 'R12')            
            n.R(1,2) = str2num(val);
        elseif strcmp(tag, 'R13')            
            n.R(1,3) = str2num(val);
        elseif strcmp(tag, 'R21')            
            n.R(2,1) = str2num(val);
        elseif strcmp(tag, 'R22')            
            n.R(2,2) = str2num(val);
        elseif strcmp(tag, 'R23')            
            n.R(2,3) = str2num(val);
        elseif strcmp(tag, 'R31')            
            n.R(3,1) = str2num(val);
        elseif strcmp(tag, 'R32')            
            n.R(3,2) = str2num(val);
        elseif strcmp(tag, 'R33')            
            n.R(3,3) = str2num(val);
        end
    end
    child = node.getFirstChild;
    b = 0;
    n.nbr = zeros(1,3)-2;
    while ~isempty(child)
        if child.getNodeType == child.ELEMENT_NODE
            if strcmp(char(child.getTagName), 'Neighbors')
                b = b+1;
                childAttribs = child.getAttributes;
                for j = 0:childAttribs.getLength-1
                    attrib = childAttribs.item(j);
                    tag = char(attrib.getName);
                    val = char(attrib.getValue);
                    if strcmp(tag, 'ID')
                        n.nbr(b) = str2num(val);
                    end
                end
            end
        end
        child = child.getNextSibling;
    end
    tracegraph(k+1) = n;
end