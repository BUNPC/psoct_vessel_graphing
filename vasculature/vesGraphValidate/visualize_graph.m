function visualize_graph(nodes, edges, title_str, node_highlight)
%visualize_graph plot graph & highlight nodes/edges in loops
% INPUTS:
%   nodes (array): node positions
%   edges (array): edges (start/end nodes)
%   title_str (string): title of chart
%   node_highlight (array): node indices to highlight

% Copy edges into standard format
s = edges(:,1); % source node
t = edges(:,2); % target node

% Create standard Matlab graph
g = graph(s, t);

%%% Plot graph
figure;
p = plot(g, 'XData', nodes(:,1), 'YData', nodes(:,2), 'ZData', nodes(:,3));
% Set nodes red
p.NodeColor = 'k';
% Set line width of edges
p.LineWidth = 2;
p.EdgeColor = [0.5 0.5 0.5];
p.MarkerSize = 3;

%%% Highlight Loops
% Determine if the graph contains cycles
[~, edgecycles] = allcycles(g);
% If so, then highlight them
if ~isempty(edgecycles)
    % Highlight edges
    for ii=1:length(edgecycles)
        highlight(p,'Edges',edgecycles{ii},'EdgeColor','r',...
                  'LineWidth',4,'NodeColor','r','MarkerSize',6)
    end
end

%%% Highlight end points and nodes in edges connecting end points
if ~isempty(node_highlight)
    % Highlight nodes
    highlight(p,node_highlight,'NodeColor','g')
end

%%% Format figure
% initialize camera view
view(3);
% Labels, title, fontsize, grid
title(title_str); xlabel('x'); ylabel('y'); zlabel('z')
set(gca, 'FontSize', 25);
grid on;
end