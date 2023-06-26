clear; close all; clc;

% load example binary skeleton image
load('skel.mat')

w = size(skel,1);
l = size(skel,2);
h = size(skel,3);

% initial step: condense, convert to voxels and back, detect cells
[A,node,link] = Skel2Graph3D(skel,0);

% total length of network
wl = sum(cellfun('length',{node.links}));

skel2 = Graph2Skel3D(node,link,w,l,h);
[A2,node2,link2] = Skel2Graph3D(skel2,0);

% calculate new total length of network
wl_new = sum(cellfun('length',{node2.links}));

% iterate the same steps until network length changed by less than 0.5%
while(wl_new~=wl)

    wl = wl_new;   
    
    skel2 = Graph2Skel3D(node2,link2,w,l,h);
    [A2,node2,link2] = Skel2Graph3D(skel2,0);
    
    wl_new = sum(cellfun('length',{node2.links}));

end

% display result
figure();
hold on;
for i=1:length(node2)
    x1 = node2(i).comx;
    y1 = node2(i).comy;
    z1 = node2(i).comz;
    
    if(node2(i).ep==1)
        ncol = 'c';
    else
        ncol = 'y';
    end
    
    for j=1:length(node2(i).links)    % draw all connections of each node
        if(node2(node2(i).conn(j)).ep==1)
            col='k'; % branches are black
        else
            col='r'; % links are red
        end
        if(node2(i).ep==1)
            col='k';
        end

        
        % draw edges as lines using voxel positions
        for k=1:length(link2(node2(i).links(j)).point)-1            
            [x3,y3,z3]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k));
            [x2,y2,z2]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k+1));
            line([y3 y2],[x3 x2],[z3 z2],'Color',col,'LineWidth',2);
        end
    end
    
    % draw all nodes as yellow circles
    plot3(y1,x1,z1,'o',...
        'MarkerFaceColor',ncol,...
        'Color','k');
end
set(gcf,'Color','white');
drawnow;
view(3);

%% Convert this graph to Matlab Graph

% The "link" struct states which nodes are connected
% s = vertcat(link2.n1);
% t = vertcat(link2.n2);
% g = graph(s,t);
g = graph(A);

% Copy node coordinates
x = vertcat(node.comx);
y = vertcat(node.comy);
z = vertcat(node.comz);

% Plot graph for comparison
figure; 
p = plot(g, 'XData',y,'YData',x,'ZData',z,'NodeColor','r');
view(3)

% Extract node locations [x,y,z]
g_x = p.XData;
g_y = p.YData;
g_z = p.ZData;

%%% Test removing nodes
% Remove node

% Plot before

% Plot after


% Calculate length of segment








