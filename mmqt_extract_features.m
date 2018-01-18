function mmqt_extract_features(trunk, figureScaling, figureHide)
% function mmqt_extract_features(fnameCZI, figureScaling, figureHide)
% function mmqt_extract_features(fnameMat, figureScaling, figureHide)
%
% Extract morphological features for microglia cells.
% This script has to be run after script mmqt_segment_image.m and mmqt_skeletonize_microglia.m
% 
% Input arguments:
%   fnameCZI      - path to CZI-file; note however, that this path-name
%                   will only be used to reconstruct the path-name of the MAT-file with the raw data (see below).
%                   In order to be found, the correspinding MAT-file should be stored in the same folder as the CZI-file.
%   fnameMat      - path to MAT-file, which contains a 4D matrix with the raw image matrix and a
%                   structure with the information about image size and scaling.
%                   This MAT-file can be created by reading CZI-files using the script: mmqt_read_convert_czi_files.m
%   figureScaling - number, indicating how big the figures should be displayed on the screen;
%                   relevant only for the figures showing orthogonal slices through the stack;
%                   a value of one means that one pixel on the screen correponds to one pixel in the image
%                   (optional, defaults to figureScaling=2.5)
%   figureHide    - logical, indicating whether to make figures visible or not
%                   If figureHide=true, figures will be made invisible and used only for saving picturs of the screenshot;
%                   (optional, defaults to figureHide=fasle)


%%
fprintf('\n\nExtracting morphological features:\n\n')

%% decompose pathname
[folder, basename, ~] = fileparts(trunk);
basename = regexprep(basename, '_stack1_raw$', ''); %--- in case the MAT-file is given, remove the suffix, which is added by mmqt_read_convert_czi_files.m
trunk = fullfile(folder, basename);

%% declare default values for unspecified variables
varNames = {'figureScaling', 'figureHide'};
varDefault = {'2.5', 'false'};
for i=1:length(varNames)
    if ~exist(varNames{i}, 'var') || isempty(eval(varNames{i}))
        eval([varNames{i} '= ' varDefault{i} ';'])
    end
end

%% Create log file
hTimer = tic;
fnameLog = [trunk, '_log3_feature_extraction.txt'];
fid = fopen(fnameLog, 'w');
fprintf(fid, '%s\n', date);
fclose(fid);

%% define how many cells to plot (i.e. plotting rule)
%--- this rule will be tested using "eval"
% plottingRule = 'mod(i,5)==0'; %--- plot every 5th cell
plottingRule = 'false'; %--- plot none of the cells
% plottingRule = 'true'; %--- plot all cells
% plottingRule = 'i==18';

%% open the watershed based parcelation
fnameAreas = [trunk,'_stack4_watershed_areas.mat'];
[hdr, MLWL] = read_3D_image_matrix(fnameAreas);
%--- remove the watershed line (ML)
WLidx = max(MLWL(:));
ML = MLWL;
ML(ML==WLidx) = 0;
MLn = WLidx-1; %--- remove one, because the maximum is the watershed-line
%--- create mask of the entire glia cells
M = MLWL>0;
%--- declare volume, to store the indices of segregated cells
MCA = zeros(size(M));

%% load the skeleton variable E,D and T (i.e. Edges, Distances, Table with node properties)
fnameMat = [trunk,'_stack4_watershed_graph.mat'];
load(fnameMat)
T.Label = (1:height(T))';
%% CREATE THE GRAPH (representing the skeleton)
G = graph(E(:,1), E(:,2), D, MLn);

%% load original segmentation
fnameSegmentation = [trunk,'_stack3_segmented.mat'];
[hdrSeg, imgSeg] = read_3D_image_matrix(fnameSegmentation);


%% load information about slices of the z-stack being included
fnameMat = [trunk,'_prepro1_spatial_corrZ_sliceOk.mat'];
temp = load(fnameMat);
sliceRange = temp.idxSliceOk;

%% extract different cell segments (whole cell, soma, nucleus)
%--- binary image of "dilated" soma
soma = zeros(size(imgSeg), 'uint8');
soma(imgSeg>=2)=1;
%--- label dileated soma
[somaL, somaN] = bwlabeln(soma);
%--- binary image of nucleus
nucleus = zeros(size(imgSeg), 'uint8');
nucleus(imgSeg==4)=1;
%--- label nucleus
[nucleusL, nucleusN] = bwlabeln(nucleus);
%---
%%
tee(fnameLog, 'Found %d soma in stack\n',somaN)
tee(fnameLog, 'Found %d nuclei in stack\n',nucleusN)

%% find for each soma the node with the largest overlap (SNV: soma-id, node-id, voxel-count)
[SNVt, ~, iSN] = unique([somaL(soma>0) ML(soma>0)],'rows');
%--- calculate number of voxel for each soma-node pair
for i=1:size(SNVt,1)
    SNVt(i,3) = sum(iSN==i);
end
%--- remove overlaps with watershed lines (which are merged to the background; label=0)
idxBackground = SNVt(:,2)==0;
SNVt(idxBackground,:) = [];
%--- indicating for each node whether it is part of a soma (overlapping with a soma)
T.inSoma = false(height(T),1);
T.inSoma(SNVt(:,2)) = true;
%--- find the node with the largest overlap
SNV = nan(somaN,3);
for i=1:somaN
    idx1 = find(SNVt(:,1)==i);
    if isempty(idx1)
        tee(fnameLog, 'WARNING: None of the nodes is overlapping with soma number %d\n', i)
        SNV(i,1) = i;
    else
        [~, idx2] = max(SNVt(idx1, 3));
        SNV(i,:) = SNVt(idx1(idx2),:);
    end
end

%% check whether multiple soma share the same seed-node, and in case merge these soma
[seedNodes, ~, idxSomaNew] = unique(SNV(:,2), 'stable');
nSeeds = length(seedNodes);
mergedSoma = false(nSeeds,1);
if nSeeds < somaN
    %--- re-number soma in order to merge soma and make numbering of soma continuous again
    overlap = zeros(nSeeds,1);
    for i=1:somaN
        overlap(idxSomaNew(i)) = overlap(idxSomaNew(i)) + SNV(i,3);
        if idxSomaNew(i) < i
            somaL(somaL==i) = idxSomaNew(i);
        end
    end
    somaN = nSeeds;
    %--- indicate which soma resulted from merging
    for i=1:nSeeds
        mergedSomaT = find(idxSomaNew==i);
        if length(mergedSomaT)>1
            mergedSoma(i) = length(mergedSomaT);
            tee(fnameLog, 'WARNING: Following soma share the same seed-node:')
            tee(fnameLog, ' %d', mergedSomaT)
            tee(fnameLog, '\n')
        end
    end
    tee(fnameLog, 'Number of merging operations (n=%d)\n', sum(mergedSoma))
    tee(fnameLog, 'Number of new soma (n=%d)\n', somaN)
    SNV = [(1:somaN)', seedNodes, overlap];
end

%% Summarize soma/seed-node relationships in a table "Cells"
Cells = array2table(SNV, 'VariableNames', {'soma','seedNode','overlap'});
Cells.mergedSoma = mergedSoma;
%--- add a variable to table T indicating whether a node is a seed (closest node to soma)
T.isSeed = false(height(T),1);
T.isSeed(seedNodes) = true;

%% Add soma specific information to the table "Cells"
S = regionprops('table', somaL,{'Centroid','Area'});
S.cogI = S.Centroid(:,[2 1 3]); %--- x/y dimension are exchanged by "regionprops"
S.cog = (S.cogI .* repmat(hdr.pixdim(2:4), height(S), 1)) - repmat(hdr.pixdim(2:4)./2, height(S), 1);
S.Centroid = [];
S.vol = S.Area .* prod(hdr.pixdim(2:4));

%--- Center of gravity (COG)
%- voxel coordinates, relative to the cropped stack
Cells.somaCogV = S.cogI;
%- voxel coordinates, relative to the original un-cropped stack (rounded)
Cells.somaCogO = round(S.cogI); %--- round to get a integer for indexing into the stack
Cells.somaCogO(:,3) = Cells.somaCogO(:,3) + sliceRange(1)-1; %--- add the sices, which were cropped due to bad quality
%- metric (microns) coordinates, relative to the cropped stack
Cells.somaCogM = S.cog;

%--- volume of soma
%- voxel count
Cells.somaVolV = S.Area;
%- metric (microns^3)
Cells.somaVol = S.vol;
Cells.branchVol = zeros(height(S),1); %--- branch volume will be calculated at the end, when the volume of the whole cell is known

%--- distance of soma COG from stack borders
%- distance from lower border (3columns for x,y,z) AND distance from upper border (3columns for x,y,z)
temp = [S.cogI - ones(height(S),3), repmat(size(ML),height(S),1) - S.cogI];
temp = temp + ones(height(S),6)./.5; %- add half a voxel, following the logic that the index gives the center of the voxel. Thus, at the border voxel the distance is still half a voxel
%- distance as voxel count
Cells.stackDistV = temp;
%- distacne in microns
Cells.stackDistM = temp .* repmat(hdr.pixdim(2:4), height(Cells),2);

%% Get number of nuclei included in each soma (SNc: soma-id, nucleus-id, voxel-count)
SNcV = unique([somaL(soma>0) nucleusL(soma>0)],'rows');
%--- remove overlap of soma with areas not being part of any nucleus
idxBackground = SNcV(:,2)==0;
SNcV(idxBackground,:) = [];
%--- get number of nuceli included in each soma 
Cells.nucleiN = zeros(somaN,1);
for i=1:somaN
    Cells.nucleiN(i) = length(SNcV(SNcV(:,1)==i, 2));
end
%---
tee(fnameLog, 'Soma with multiple nuclei inside (n=%d)\n', sum(Cells.nucleiN>1))

%% Find for each "node" in skeleton (i.e. ML) the closest "Seed-node" and "Seed-ID"
distNodes2Seeds = distances(G,seedNodes);
[~, seedIds] = min(distNodes2Seeds,[],1);
%--- find nodes not connected to any seed-node (the shortest path to all see-nodes is in this case Inf)
idNotInCell = find(all(isinf(distNodes2Seeds)));
tee(fnameLog, 'Nodes not connected to any seed-node (n=%d)\n', length(idNotInCell))
seedIds(idNotInCell) = 0; %--- set seed-node-id equal "0" for nodes not connected to any seed node
%---
T.idCell = seedIds';

%% find bridges between cells
Eid = nan(size(E,1),1);
Ecolor = nan(size(E,1),1);
BridgeCells = zeros(0,2); %--- each row contains the IDs of the cells connected by a bridge
BridgeNodes = zeros(0,2); %--- each row contains the IDs of the nodes connected by a bridge
for i=1:size(E,1)
    %--- get for current edge the seed IDs of the two connected nodes
    cellIDT = seedIds(E(i,:));
    %--- check whether these IDs belong to the same cell
    if cellIDT(1)==cellIDT(2)
        Eid(i) = cellIDT(1);
        if cellIDT(1)==0 %--- not belonging to any cell, as far as defined by a soma
            Ecolor(i) = 1;
        else
            Ecolor(i) = 2;
        end
    else %--- this means that the edge connects two different cells
        Eid(i) = 0;
        Ecolor(i) = 3;
        BridgeCells = [BridgeCells; cellIDT];
        BridgeNodes = [BridgeNodes; E(i,:)];
    end
end
%---
T.isBridgeNode = false(height(T),1);
T.isBridgeNode(BridgeNodes(:)) = true;

%% find groups of interconnected cells
%--- remove duplicates
[BridgeCellsU, ~, idxBridgeCells] = unique(sort(BridgeCells,2),'rows');

%--- use "graph" to find groups of connected cells
LC = graph(BridgeCellsU(:,1),BridgeCellsU(:,2), [], nSeeds);
Cells.group = conncomp(LC)'; %, 'OutputForm', 'cell');
%--- for each seed collect the degree of conncetedness in its group
Cells.groupDegree = zeros(nSeeds,1);
Cells.groupDegreeMultiple = zeros(nSeeds,1);
for i=1:size(BridgeCellsU,1)
    Cells.groupDegree(BridgeCellsU(i,:))         = Cells.groupDegree(BridgeCellsU(i,:))         + 1;
    Cells.groupDegreeMultiple(BridgeCellsU(i,:)) = Cells.groupDegreeMultiple(BridgeCellsU(i,:)) + sum(idxBridgeCells==i);
end

%% For each seed get the shortest path to the closest connected seed
distSeed2Seed = distances(G,seedNodes, seedNodes);
distSeed2Seed(eye(nSeeds)>0) = NaN; %--- diagonal indicates the distance of a cell to itself, which is Zero
% distSeed2Seed(isinf(distSeed2Seed)) = NaN; %--- "Inf" could be substituted with "NaN", given that an "Inf" indicate that two cells are not connected at all
Cells.groupShortestDist = min(distSeed2Seed, [], 'omitnan')';


%% declare table columns
Cells.stackTouch = zeros(height(Cells),6);
% Seeds.stackDist = zeros(height(Seeds),6);
Cells.cellVol = zeros(height(Cells),1);
Cells.cellArea = zeros(height(Cells),1);
Cells.somaNodesVol = zeros(height(Cells),1);
Cells.branchNodesVol = zeros(height(Cells),1);
Cells.sphericity = zeros(height(Cells),1);
Cells.circularity = zeros(height(Cells),1);
Cells.solidity = zeros(height(Cells),1);
Cells.nodesN = zeros(height(Cells),1);
Cells.nodesBranching = zeros(height(Cells),1);
Cells.nodesEnding = zeros(height(Cells),1);
Cells.nodesBranchingRatio = zeros(height(Cells),1);
Cells.nodesEndingRatio = zeros(height(Cells),1);
%---
Cells.degreeSeed = zeros(height(Cells),1);
Cells.closenessSeed = zeros(height(Cells),1);
Cells.betweennessSeed = zeros(height(Cells),1);
Cells.volumeSeed = zeros(height(Cells),1);
%---
Cells.degree = zeros(height(Cells),5);
Cells.closeness = zeros(height(Cells),5);
Cells.betweenness = zeros(height(Cells),5);
Cells.volume = zeros(height(Cells),5);
%---
Cells.branchN = zeros(height(Cells),1);
Cells.branchNodesN = zeros(height(Cells),1);
Cells.branchEndNodesN = zeros(height(Cells),1);
Cells.branchEndNodesPerBranch = zeros(height(Cells),1); %--- corresponds to ramification index
Cells.branchSegmentsN = zeros(height(Cells),1);
Cells.branchCyclesN = zeros(height(Cells),1);
Cells.branchNodesPerSegment = zeros(height(Cells),1);
Cells.branchNodesPerBranch = zeros(height(Cells),1); 
Cells.branchSegmentsPerBranch = zeros(height(Cells),1);
Cells.branchLengthSkel = zeros(height(Cells),5);
Cells.branchLengthAir = zeros(height(Cells),5);
Cells.branchLengthRatio = zeros(height(Cells),5);


%--- declare cell array to store the graph for each cell
cellGraphs = cell(height(Cells),3);


%% segregate cells
for i=1:nSeeds
    tic
    close all
    %%
    tee(fnameLog, '\nAnalysing cell %d out of %d cells\n', i, nSeeds)
    %% extract subgraph of nodes with closest connection to seed
    nodeIDs = find(seedIds==i);
    TC = T(nodeIDs,:);
    GC = subgraph(G, nodeIDs);
    GC.Nodes.Name =cellfun(@num2str,table2cell(T(nodeIDs,'Label')),'Uniform',false);
    GC.Nodes.Index = T{nodeIDs,'Label'};
    
    %% plot the graph
    if eval(plottingRule)
        hf = figure;
        hf.Position = [1292 120 1265 1225];
        hf.Color = 'white';
        ha = axes(); hold on
        % hg =plot(GS);
        hg =plot(GC,'XData',T.cog(nodeIDs,1),'YData',T.cog(nodeIDs,2),'ZData',T.cog(nodeIDs,3));
        %hg =plot(G,'XData',T.depth(:,1),'YData',T.depth(:,2),'ZData',T.depth(:,3));
        adjust_viewing_options(ha, hdr, T(nodeIDs,:))
        adjust_graph_properties(hg, T(nodeIDs,:))
    end
    %% Node count
    Cells.nodesN(i) = height(TC);
    Cells.nodesBranching(i) = sum(TC.degree>2);
    Cells.nodesEnding(i) = sum(TC.degree==1);
    Cells.nodesBranchingRatio(i) = sum(TC.degree>2)/height(TC);
    Cells.nodesEndingRatio(i) = sum(TC.degree==1)/height(TC);
    %% Node Degree
    Cells.degreeSeed(i) = TC.degree(TC.isSeed);
    Cells.degree(i,1:5) = prctile(TC.degree,[0,25,50,75,100]);
    %% Node Closeness
    TC.closeness = centrality(GC,'closeness');
    Cells.closenessSeed(i) = TC.closeness(TC.isSeed);
    Cells.closeness(i,1:5) = prctile(TC.closeness,[0,25,50,75,100]);
    %% Node Betweenness
    TC.betweenness = centrality(GC,'betweenness');
    Cells.betweennessSeed(i) = TC.betweenness(TC.isSeed);
    Cells.betweenness(i,1:5) = prctile(TC.betweenness,[0,25,50,75,100]);
    %% Node Volume
    Cells.volumeSeed(i) = TC.vol(TC.isSeed);
    Cells.volume(i,1:5) = prctile(TC.vol,[0,25,50,75,100]);
    %% Volume of soma and branches (based on the volume of the nodes overlapping with the soma)
    Cells.somaNodesVol(i) = sum(TC.vol(TC.inSoma));
    Cells.branchNodesVol(i) = sum(TC.vol(~TC.inSoma));

    %% Create "digraph" starting at the seed-node
    seedNode = find(TC.isSeed);
    GCD = shortestpathtree(GC, seedNode); %, 'Method', 'unweighted')
    %--- Find lost edges (edges are lost, if there are cycles in the graph)
    if size(GCD.Edges,1)<size(GC.Edges,1)
        %--- create copy of GCD whith edges building cycles
        GCDC = GCD;
        %--- convert edge.EndNode names in GC and GCD (i.e. indices in G) from char to numbers
        GCedges = cellfun(@str2double, GC.Edges.EndNodes);
        GCDedges = cellfun(@str2double, GCD.Edges.EndNodes);
        %--- find lost edges
        idxEdgeMissing = [];
        for j=1:size(GCedges,1)
            if ~any(any(GCDedges==GCedges(j,1),2) & any(GCDedges==GCedges(j,2),2))
                idxEdgeMissing = [idxEdgeMissing, j];
            end
        end
        %--- reinsert lost edges
        for j=1:length(idxEdgeMissing)
            d = distances(GC, findnode(GC,seedNode), GC.Edges.EndNodes(idxEdgeMissing(j),:));
            if d(1)<=d(2)
                GCDC = addedge(GCDC, GC.Edges.EndNodes(idxEdgeMissing(j),1), GC.Edges.EndNodes(idxEdgeMissing(j),2), GC.Edges.Weight(idxEdgeMissing(j)));
            else
                GCDC = addedge(GCDC, GC.Edges.EndNodes(idxEdgeMissing(j),2), GC.Edges.EndNodes(idxEdgeMissing(j),1), GC.Edges.Weight(idxEdgeMissing(j)));
            end
        end
    else
        %--- clear the variable to avoid using the graph from the previous cell
        clearvars GCDC
    end
    %--- find indices of edges.EndNodes
    Edig = reshape(findnode(GCD, GCD.Edges.EndNodes), [],2);
    
    %% find major branches, defining them as subgraphs starting from a node outside the soma and a length > 2 microns
    %--- find all first non-soma nodes as potential branch start-nodes
    branchStart = unique(Edig(TC.inSoma(Edig(:,1)) & ~TC.inSoma(Edig(:,2)) , 2));
    %--- find all dead-end nodes as potential branch end-nodes
    branchEnd = find(TC.degree==1 & ~TC.isSeed);
    %--- calculate the distance between each potential start and end point
    branchDist = distances(GCD, branchStart, branchEnd);
    branchDist(isinf(branchDist)) = nan;
    %--- find for each potential start node the end-node with the longest distance
    [branchDist, idxTemp] = max(branchDist, [],2, 'omitnan');
    branchEndT = branchEnd(idxTemp);
    branches = [branchStart, branchEndT, branchDist];
    %--- exclude branches which are shorter than 2 microns
    branches = branches(branchDist > 2 , :);

    %---
    Cells.branchN(i) = size(branches,1);
    %---
    if Cells.branchN(i)>0
        Cells.branchLengthSkel(i,:) = prctile(branches(:,3),[0,25,50,75,100]);
        branchLengthAir = sum((T.cog(GCD.Nodes.Index(branches(:,1)),:) - T.cog(GCD.Nodes.Index(branches(:,2)),:)).^2, 2).^.5;
        Cells.branchLengthAir(i,:) = prctile(branchLengthAir,[0,25,50,75,100]);
        %---
        branchLengthRatio = branches(:,3)./branchLengthAir;
        Cells.branchLengthRatio(i,:) = prctile(branchLengthRatio,[0,25,50,75,100]);
        %--- count all end nodes for all branches
        branchDist = distances(GCD, branches(:,1), branchEnd);
        Cells.branchEndNodesN(i) = sum(~all(isinf(branchDist)));
        Cells.branchEndNodesPerBranch(i) = sum(~all(isinf(branchDist)))/size(branches,1); %--- ramification index
    else
        Cells.branchEndNodesPerBranch(i) = NaN; %--- if there are no branches, the ratio is undefined (devision by zero in Matlab gives Inf instead, which we don't want)
    end
    %% analyze subgraph for each major branch
    branchesEdges = cell(0);
    for j=1:size(branches,1)
        %%
        GCDB = shortestpathtree(GCD, branches(j,1));
        %--- remove non-connected nodes
        idxNonConnected = find( (indegree(GCDB) + outdegree(GCDB)) == 0);
        GCDB = rmnode(GCDB,idxNonConnected);
        
        %--- Given that the shortest-path-tree looses edges for cyclic
        %- graphs, extract the subgraph for the branch from the di-graph including cycles (i.e. GCDC), using all nodes of the
        %- shortest-path-tree
        if exist('GCDC', 'var')
            GCDB = subgraph(GCDC, GCDB.Nodes.Name);
        end

        %--- extract edges from digraph (by name) 
        %- for removal of edges leading back into the soma (might happen for cyclic graphs)
        branchEdgesNames = GCDB.Edges.EndNodes;
        branchesEdges{j} = branchEdgesNames;
        idx1 = cellfun(@str2num, reshape(branchEdgesNames,[],1));
        idx2 = T.inSoma(idx1);
        if any(idx2)
            nodeNames = arrayfun(@num2str,idx1(idx2),'UniformOutput', false);
            GCDB = rmnode(GCDB, nodeNames);
        end
        %--- extract edges from digraph (by name) 
        %- for highlighting branches in the graph plot
        if eval(plottingRule)
            branchEdgesNames = GCDB.Edges.EndNodes;
            highlight(hg, branchEdgesNames(:,1), branchEdgesNames(:,2), 'EdgeColor', 'c')
        end
        
        %--- determine number of branches going in and out of each node (bio)
        bio = [indegree(GCDB), outdegree(GCDB)];
        % bio(:,3) = bio(:,1) + bio(:,2);

        %--- find number of edges leading away (i.e. outdegree or bio(:,2)) from a
        %- non-reachable node or from nodes with more than one in- or out-going edge.
        %- the sum of these edges corresponds to the number of branch segments
        Cells.branchSegmentsN(i) = Cells.branchSegmentsN(i) + sum(bio(bio(:,1)~=1 | bio(:,2)>1, 2));
        
        %--- find nodes with indegree > 1; the indegree-1 of these nodes corresponds to the number of cycles in the graph
        idxConverge = find(bio(:, 1)>1);
        Cells.branchCyclesN(i) = Cells.branchCyclesN(i) + sum(bio(idxConverge, 1)) - length(idxConverge);
        
        %--- if a node has exactly indegree==2 and outdegree==0, then the two incoming branch-segments shold be counted as one only
        %--- => subtract count of these nodes from branchSegmentsN
        Cells.branchSegmentsN(i) = Cells.branchSegmentsN(i) - sum(bio(idxConverge, 1)==2 & bio(idxConverge, 2)==0);

        %--- number of nodes in the branch
        Cells.branchNodesN(i) = Cells.branchNodesN(i) + height(GCDB.Nodes);
        
    end
    %---
    Cells.branchNodesPerSegment(i) = Cells.branchNodesN(i)/Cells.branchSegmentsN(i);
    Cells.branchNodesPerBranch(i) = Cells.branchNodesN(i)/Cells.branchN(i);
    Cells.branchSegmentsPerBranch(i) = Cells.branchSegmentsN(i)/Cells.branchN(i);
    %%
    cellGraphs{i,1} = GC;
    cellGraphs{i,2} = T(nodeIDs,:);
    cellGraphs{i,3} = branchesEdges;
    
    %% create a mask containing only the cell
    tic
    %--- remove all nodes not being part of the cell
    MX = MLWL;
    MX(~ismember(MLWL,nodeIDs)) = 0;
    %--- dilate the areas and mask with M, in order to reinsert the watershed-line voxel
    MX = imdilate(MX>0, ones([3,3,3])) & M;
    %--- check whether the dilation swept into non-included nodes and remove theses voxel
    voxelNotInCellN = sum(~ismember(MLWL(MX), [nodeIDs WLidx]));
    if voxelNotInCellN>0
        tee(fnameLog, 'Dilation swept into areas of nodes, not belonging to the current cell (#voxel=%d)\n', voxelNotInCellN)
        MX(~ismember(MLWL, [nodeIDs WLidx])) = false;
    end
    toc
    %--- check if the mask contains more than one connected area (and keep only the biggest area) NOTE: This should hopefully not happen!
    [MXL, N] = bwlabeln(MX,26);
    if N>1
        MXT = regionprops('table', MXL, 'Area');
        tee(fnameLog, 'Dilation flew into %d non-connected area(s), which will be removed\n', N-1)
        tee(fnameLog, 'Size of areas:\n')
        tee(fnameLog, '   %d\n', MXT.Area)
        [~, idx] = max(MXT.Area);
        MX(MXL~=idx) = 0;
    end
    %--- determine the range of voxel containing the whole cell, for cropping the volume
    [X, Y, Z] = ind2sub(size(MX), find(MX));
    cropRange = [min(X) max(X); min(Y) max(Y); min(Z) max(Z)];
    
    %% update mask containing cell-IDs for all voxel belonging to that cell
    MCA(MX) = i;
    
    
    %% Check if the cell touches the border of the image stack and add a voxel to the crop volume otherwise
    cellTouches = [cropRange(:,1)==[1; 1; 1], cropRange(:,2)==size(ML)'];
    cropRange = cropRange + ~cellTouches .* [-1,1; -1,1; -1,1];
    
    %% Coordinates should follow the order (i.e. x-lower-margin, y-lower-margin, z-lower-margin, x-upper-margin, y-upper-margin, z-upper-margin)
    Cells{i,'stackTouch'} = reshape(cellTouches,1,6);
    %% crop - general mask
    MC = MX(cropRange(1,1):cropRange(1,2), cropRange(2,1):cropRange(2,2), cropRange(3,1):cropRange(3,2));
       
    %% determine volume and surface area of cell
    [~, voxT, ~, volT, areaT] = find_individual_clusters_find_surface_anisotrop(MC,0.5, hdr.pixdim(2:4));
    %--- given the low connectivity (conn=6) used in this function, more than one cluster can result 
    tee(fnameLog, '- detected %d cluster with (conn=6)\n', length(voxT))
    %tee(fnameLog, '  -> %d/%d cluster have 5 or less voxel\n', length(voxT), sum(voxT<=5))
    tee(fnameLog, '  voxel count:')
    tee(fnameLog, ' %d', voxT); tee(fnameLog, '\n');
    %---
    areaT = sum(areaT);
    volT = sum(volT);
    Cells.cellArea(i) = areaT;
    Cells.cellVol(i) = volT;  
    
    %% calculate sphericity index
    Cells.sphericity(i) = pi^(1/3) * (6*volT).^(2/3) ./ areaT; %--- Wikipedia article on "Sphericity" (Defined by Wadell in 1935)
    
    %% project in z-direction and calculate circularity and solidity
    MCP = any(MC,3);
    % MCP = imfill(MCP, 'holes'); %--- fill holes to avoid discontiguous regions?
    temp = regionprops(MCP,{'Area','Perimeter','Solidity'});
    Cells.circularity(i) = 4*pi*temp.Area / temp.Perimeter^2;
    Cells.solidity(i) = temp.Solidity;
    
    %% iso-surface
    if eval(plottingRule)
        add_isosurface(MC, hdr, cropRange(:,1)-1, cellTouches);
        drawnow
    end
    %% crop labeled volume for video
    MLWLC = MLWL(cropRange(1,1):cropRange(1,2), cropRange(2,1):cropRange(2,2), cropRange(3,1):cropRange(3,2));
    MOther = MLWLC;
    MOther(MC>0) = 0;
    MOther(MOther>0) = 1;
    MLWLC(~MC) = 0;
        
    %%
    if eval(plottingRule)
        make_video_single_cell(hf, sprintf('%s_cell%03d_branches_test_end_nodes', trunk, i), MLWLC, hdr, T(nodeIDs,:), cropRange(:,1)-1, MOther);
    end
    %%
    display_time_delay()
end

%% calculate volume of branches
Cells.branchVol = Cells.cellVol - Cells.somaVol;



%% SAVING RESULTS
%% save the table
writetable(Cells, [trunk '_features_of_cells.txt'], 'Delimiter', '\t');

%% final result of segregation
% cmap = [0 0 0; jet(max(MCA(:)))];
cmap = [0.5 0.5 0.5; jet(max(MCA(:)))];
hdr.cmap = cmap;
%--- show final result
h = show_cells_3D(MCA,hdr,[],[],figureScaling, [], figureHide);

%--- save picture
F = getframe(h);
imwrite(F.cdata, [trunk,'_ortho5_segregated_cells.png']);


%% save image containing segregated cells
fnameOut = [trunk, '_stack5_segregated_cells.mat'];
hdr = change_hdr_dim(hdr, size(MCA));
if height(Cells) <= intmax('uint8')
    hdr = change_hdr_datatype_info(hdr, 2);
    write_3d_image_matrix(hdr, uint8(MCA), fnameOut)
else
    hdr = change_hdr_datatype_info(hdr, 4);
    write_3d_image_matrix(hdr, int16(MCA), fnameOut)
end
    
%% save graphs for each single cell
fnameOut = [trunk, '_stack5_segregated_graphs.mat'];
save(fnameOut, 'cellGraphs');

%% save the (extended) skeleton variables E,D and T (i.e. Edges, Distances, Table with node properties)

fnameMat = [trunk,'_stack4_watershed_graph_extended.mat'];
save(fnameMat, 'E', 'D', 'T')


%%
tee(fnameLog, 'Feature extraction completed!\n')
tee(fnameLog, display_time_delay(hTimer))
fprintf('\n\n')
