%% Plot graph with lines
function scatter_graph(edges, nodes, tstring)
% Scatter plot from edges and nodes
figure;
% Scatter plot of nodes
scatter3(nodes(:,1), nodes(:,2), nodes(:,3), '.','b');
% Plot edges
for ii=1:length(edges)
    n1 = edges(ii,1);
    n2 = edges(ii,2);
    x = [nodes(n1,1), nodes(n2,1)];
    y = [nodes(n1,2), nodes(n2,2)];
    z = [nodes(n1,3), nodes(n2,3)];
    line(x, y, z, 'Color', 'red');
    title(tstring);
end
end