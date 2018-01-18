function adjust_graph_properties(hg, T)
%% adjust graph properties
%--- define node size (i.e. marker size)
Gdeg = T.degree;
MarkerSize = Gdeg;
MarkerSize(Gdeg==0) = 10; %--- unconnected nodes
MarkerSize(Gdeg==2) = 5; %0.001; %--- non-branching nodes
MarkerSize(Gdeg==1) = 7;    %--- start/end nodes
MarkerSize(Gdeg>2) = 7;     %--- branching points
if any(strcmp('isSeed', T.Properties.VariableNames))
    MarkerSize(T.isSeed) = 15;     %--- nucleus
    MarkerSize(T.isBridgeNode) = 12;
end
% MarkerSize(idCell==0) = 0.001;
% hg.Marker = 'o';
% hg.MarkerSize = 5;
% hg.MarkerSize = 2;
hg.MarkerSize = MarkerSize;

%--- define node colors
% G.Nodes.NodeColors = distNucleusColor';
% G.Nodes.NodeColors = idCell';
% hg.NodeCData = G.Nodes.NodeColors;
cmap = repmat([0.6 0.4 0.4], height(T), 1); %--- default color brown
cmap(T.degree==1, :) = repmat([1 0 1], sum(T.degree==1), 1); %--- magenta for end points
cmap(T.degree>2, :) = repmat([0 0.8 0.3], sum(T.degree>2), 1); %--- green for branching points
if any(strcmp('inSoma', T.Properties.VariableNames))
    cmap(T.inSoma, :) = repmat([1 0.549 0], sum(T.inSoma), 1); %--- red for seed
end
if any(strcmp('isSeed', T.Properties.VariableNames))
    cmap(T.isSeed, :) = repmat([1 0 0], sum(T.isSeed), 1); %--- red for seed
end
if any(strcmp('isBridgeNode', T.Properties.VariableNames))
    cmap(T.isBridgeNode, :) = repmat([0 0 1], sum(T.isBridgeNode), 1); %--- blue for Bridges to other cells
end
hg.NodeColor = cmap;
% hg.NodeColor = [0 0.8 0.3];

%--- color map
colormap(jet)
% colorbar

%--- define edge properties
hg.LineWidth = 2;
% hg.LineWidth = 1;
% hg.LineWidth = ((Ecolor>2)*4+1);
hg.EdgeAlpha = 1;
% cmap = [0.7 0.7 0.7; 0.4 0.4 1; 1 0 1];
% cmap = cmap(Ecolor,:);
% hg.EdgeColor = cmap;
% hg.EdgeColor = [0.8 0 0.2];
hg.EdgeColor = [0.6 0.4 0.4];
