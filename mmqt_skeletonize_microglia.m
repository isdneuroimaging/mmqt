function mmqt_skeletonize_microglia(trunk)
% function mmqt_skeletonize_microglia(fnameCZI)
% function mmqt_skeletonize_microglia(fnameMat)
%
% Create skeleton of microglia cells. 
% This script has to be run after script mmqt_segment_image.m
% 
% Input arguments:
%   fnameCZI      - path to CZI-file; note however, that this path-name
%                   will only be used to reconstruct the path-name of the MAT-file with the raw data (see below).
%                   In order to be found, the correspinding MAT-file should be stored in the same folder as the CZI-file.
%   fnameMat      - path to MAT-file, which contains a 4D matrix with the raw image matrix and a
%                   structure with the information about image size and scaling.
%                   This MAT-file can be created by reading CZI-files using the script: mmqt_read_convert_czi_files.m



%%
fprintf('\n\nCreating microglia skeleton:\n\n')

%% decompose pathname
[folder, basename, ~] = fileparts(trunk);
basename = regexprep(basename, '_stack1_raw$', ''); %--- in case the MAT-file is given, remove the suffix, which is added by mmqt_read_convert_czi_files.m
trunk = fullfile(folder, basename);


%% Check existence of segmentation file
fnameSegmentation = [trunk,'_stack3_segmented.mat'];
if ~exist(fnameSegmentation,'file')
    fprintf('ERROR: The segmentation file does not exist:\n %s\n', fnameSegmentation)
    fprintf('-> Aborting skeletonization\n')
    return
end

hTimer = tic;

%% open/create file for logging info
fnameLog = [trunk '_log2_skeletonization.txt'];
fidLog = fopen(fnameLog,'w');
fprintf(fidLog, '%s\n', date);

%% read the image with the segmentation
tee(fidLog, 'Reading segmentation from file... ')
tic
[hdr, img] = read_3D_image_matrix(fnameSegmentation);
tee(fidLog, display_time_delay)

%% create mask for the whole cells
M = img>0;

%% distance transformation, as negative value
MD = bwdist(~M);
MI = -MD;

%% watershed
tee(fidLog, 'Watershed segmentation... ')
tic
conn= 18;
ML = watershed(MI, conn);
tee(fidLog, display_time_delay)
%% remove the background from the cluster
ML(~M) = 0;
MLn = double(max(ML(:)));

%%
%% get watershed lines (remaining voxel, not part of any watershed-area)
MS = M;
MS(ML>0) = 0;

%% check for all voxel on the watershed-line, which cluster they are bordering with
tee(fidLog, 'Determining neighborhood tuples of watershed-line voxel... ')
tic
%--- get indices and subscripts of each watershed-line voxel
Uind = find(MS==1);
[x,y,z] = ind2sub(size(MS), Uind);
%--- define a box around each voxel
x = [x-1, x+1];
y = [y-1, y+1];
z = [z-1, z+1];
%--- crop the box at the lower border of the matrix
x(x(:,1)<1, 1) = 1;
y(y(:,1)<1, 1) = 1;
z(z(:,1)<1, 1) = 1;
%--- crop the box at the upper border of the matrix
x(x(:,2)>size(MS,1), 2) = size(MS,1);
y(y(:,2)>size(MS,2), 2) = size(MS,2);
z(z(:,2)>size(MS,3), 2) = size(MS,3);
%--- get the neighbors for each voxel on the watershed line
U = cell(length(Uind),1);
for i=1:length(Uind)
    temp = unique(ML(x(i,1):x(i,2), y(i,1):y(i,2), z(i,1):z(i,2)));
    U{i} = temp(2:end); %--- remove the "0" which is always present as the value of the watershed line voxel itself
end
tee(fidLog, display_time_delay)


%% merge areas without contact to the background, to areas in their neighboring
%% check which area has no or very little contact to the background and merge these areas with its neighbors
tee(fidLog, 'Checking contact with background for each area... ')
tic
Mtemp = single(ML) + (single(MS)*(MLn+1));
[internal, background] = find_clusters_neighbors( Mtemp , hdr.pixdim(2:4));
background(end) =[]; %--- remove the watershed-line area which had the largest label
internal=internal(1:end-1, end);
%--- consider areas without contact to the surface (i.e. background) AND areas with little contact to the surface
idxLogical= background<0.9 & background./internal<1/3; %--- contact with surface <0.9 micron & surface/internal<1/3
idxWomb = find(idxLogical);
idxSurf = find(~idxLogical);

tee(fidLog, display_time_delay)

%% reassign labels for areas without contact to background
if ~isempty(idxWomb)
    tee(fidLog, '-> There are %d areas without or with little contact to the background\n', length(idxWomb))
    fprintf(fidLog,' %d', idxWomb);
    fprintf(fidLog,'\n');
    tee(fidLog, 'Merging these areas to their neighbors... ')
    tic
    %% find neighborhood matrix (i.e. find for each pair of areas the number of bridging watershed-line voxel)
    neighbors = zeros(MLn, MLn);
    for i=1:size(U,1)
        %--- increment the neighborhood matrix if two areas touch the same watershe-line voxel
        UT = U{i};
        for j=1:length(UT)-1
            for k=j+1:length(UT)
                neighbors(UT(j),UT(k)) = neighbors(UT(j),UT(k)) + 1;
            end
        end
    end
    %--- make neighborhood matrix symetric
    neighbors = neighbors + neighbors';
    %% for areas without contact to the background, find the neighboring area with the largest surface of contact
    labelOldNew = nan(length(idxWomb),3);
    for i = 1:length(idxWomb)
        neT = neighbors(idxWomb(i),:);
        [maxt, idx2] = max(neT(idxSurf));
        if maxt>0 %--- new area has contact to surface
            labelOldNew(i,:) = [idxWomb(i), idxSurf(idx2), 1];
        else
            [maxt, idx2] = max(neT(idxWomb));
            if maxt>0 %--- new area is itself in the womb
                labelOldNew(i,:) = [idxWomb(i), idxWomb(idx2), 0];
            end
        end
    end    
    
    %% make sure that the reassignment to non-surface cluster happens before these cluster are changed
    if ~all(labelOldNew(:,3))
        tee(fidLog, 'WARNING: there are non-surface cluster which will be merged to another non-surface cluster.\n')
        labelOldNew = sortrows(labelOldNew, 3); %--- make sure that the reassignment to non-surface cluster happens befor ethese cluster are changed
    end
    
    
    %% change lable in variable U
    nChange = size(labelOldNew,1);
    for i=1:length(Uind)
        UT = U{i};
        for j=1:nChange
            UT(UT==labelOldNew(j,1)) = labelOldNew(j,2);
        end
        U{i} = unique(UT);
    end
    
    %% get length of neighborhood tuples in U
    Ul = cellfun(@length,U);
    
    %% insert label in ML, for all watershed line voxel who are left with a single neighboring area
    ML(Uind(Ul==1)) = [U{Ul==1}];
    %%
    %% remove these watershed-line voxel with a single neighboring area from variables MS, U, Uind and Ul
    MS(ML>0) = 0;
    U(Ul==1,:) = [];
    Uind(Ul==1) = [];
    Ul(Ul==1) = [];
    
    %% change lable in variable ML
    for j=1:nChange
        ML(ML==labelOldNew(j,1)) = labelOldNew(j,2);
    end
    
    %%
    tee(fidLog, display_time_delay)
else
    %% get length of neighborhood tuples in U
    Ul = cellfun(@length,U);
end

%% Check and remove neighborhood tuples for voxel not touching any area or touching only one area (both cases should be very rare)
idx = find(Ul==0);
if ~isempty(idx)
    tee(fidLog, '- there are %d surface voxel not touching any area\n', length(idx))
    U(idx) = [];
    Uind(idx) = [];
    Ul(idx) = [];
end
idx = find(Ul==1);
if ~isempty(idx)
    tee(fidLog, '- there are %d surface voxel touching only one area\n', length(idx))
    U(idx) = [];
    Uind(idx) = [];
    Ul(idx) = [];
end
%--- sort U by length of neighborhood
[Ul, I] = sort(Ul);
U = U(I);
Uind = Uind(I);

%% Find unique neighborhood tuples (variable UN)
tee(fidLog, 'Find unique neighborhood tuples... ')
tic
%--- convert cell of number vectors of different length to zero-padded matrix 
%--- this takes time, but is needed for the "unique" function, to find unique neighborhood tuples
lvn = repmat({max(Ul)},length(U),1);
padZero = @(v, lvn) padarray(v, lvn-length(v),'post');
UZ = cellfun(padZero, U, lvn, 'UniformOutput', false);
UZ = cell2mat(UZ')';
%--- find unique 
[~, iUnique] = unique(UZ,'rows', 'stable');
UN = U(iUnique);
tee(fidLog, display_time_delay)


%%  Check for each neighborhood tuple, whether it is included in a longer neighborhood tuple 
%--- and exclude the shorter tuple in this case
%--- Given that non-unique tuples were allready removed, tuples have to be
%--- tested only against longer tuples
tee(fidLog, 'Check for each neighborhood tuple, whether it is included in a longer tuple... ')
tic
%--- get length of neighborhood tuples in U
Ul = cellfun(@length,UN);
startIdx = [1; find(diff(Ul))+1];
%---
include = true(size(UN,1),1);
for iC = 1:length(startIdx)-1
    for i=startIdx(iC):startIdx(iC+1)-1
        for j=startIdx(iC+1):size(UN,1)
            if all(ismember(UN{i}(:), UN{j}(:)))
                include(i) = false;
                break
            end
        end
    end
end
UN = UN(include);
tee(fidLog, display_time_delay)


%% find interlaced tuples in UN, and find nodes overlapping between interlaced tuples
%--- interlaced neighborhoods shell be neighborhoods which share more at
%--- least two areas. This has to be tested only for tuples of length 3 and
%--- higher, given that tuples of length 2 have been already included, if they
%--- were included in a longer one
tee(fidLog, 'Find interlaced neighborhoods and record the overlapping areas for merging... ')
tic
UNl = cellfun(@length,UN);
startIdx = find(UNl==3,1,'first');
interlaced = false(length(UN), length(UN));
merge = cell(0);
k=0;
for i=startIdx:length(UN)
    for j=i+1:length(UN)
        overlap = ismember(UN{i}(:), UN{j}(:));
        if sum(overlap) > 1
            k=k+1;
            interlaced(i,j) = true;
            merge{k} = UN{i}(overlap);
        end
    end
end
tee(fidLog, display_time_delay)


%%
tee(fidLog, 'Prepare reassignment table for merging... ')
tic
%--- check whether nodes to be merged overlap and merge these groups
for i=1:length(merge)
    for j=i+1:length(merge)
        if any(ismember(merge{i}(:), merge{j}(:)))
            merge{j} = unique(vertcat(merge{[i,j]}));
            merge{i} = [];
            break
        end
    end
end
%--- remove cells from merge which ar now empty
lm = cellfun(@length, merge);
merge(lm==0) = [];
lm(lm==0) = [];

%--- create assignment table "labelOldNew" for merging operation
labelOldNew = nan(sum(lm)-length(lm),1);
k=0;
for i=1:length(merge)
    ll = length(merge{i})-1;
    labelOldNew(k+1:k+ll,1) = merge{i}(2:end);
    labelOldNew(k+1:k+ll,2) = merge{i}(1);
    k = k+ll;
end
tee(fidLog, display_time_delay)

%% merge by changing the area labels
tee(fidLog, 'Merge by changing the area labels... ')
tic
%--- change lable in variable U
nChange = size(labelOldNew,1);
for i=1:length(U)
    UT = U{i};
    for j=1:nChange
        UT(UT==labelOldNew(j,1)) = labelOldNew(j,2);
    end
    U{i} = unique(UT);
end

%--- change lable in variable UN
nChange = size(labelOldNew,1);
for i=1:length(UN)
    UT = UN{i};
    for j=1:nChange
        UT(UT==labelOldNew(j,1)) = labelOldNew(j,2);
    end
    UN{i} = unique(UT);
end

%--- get length of neighborhood tuples in U and UN
Ul = cellfun(@length,U);
UNl = cellfun(@length,UN);

%--- give watershed-line voxel between merged areas an area-label (i.e. watershed-line voxel who are left with a single neighboring area)
ML(Uind(Ul==1)) = [U{Ul==1}];

%--- remove these watershed-line voxel with a single neighboring area from variables MS, U, Uind, Ul and in UN, UNl
MS(Uind(Ul==1)) = 0;
%---
U(Ul==1,:) = [];
Uind(Ul==1,:) = [];
Ul(Ul==1) = [];
%---
UN(UNl==1,:) = [];
UNl(UNl==1) = [];


%--- change lable in variable ML
for j=1:nChange
    ML(ML==labelOldNew(j,1)) = labelOldNew(j,2);
end
tee(fidLog, display_time_delay)
    
%% define neighborhood matrix (i.e. find for each pair of areas the number of bridging watershed-line voxel)
tee(fidLog, 'Creating neighborhood matrix for all areas... ')
tic
neighbors = zeros(MLn, MLn);
for i=1:length(U)
    UT = U{i};
    for j=1:length(UT)-1
        for k=j+1:length(UT)
            neighbors(UT(j),UT(k)) = neighbors(UT(j),UT(k)) + 1;
        end
    end
end
%--- make the matrix symetric
neighbors = neighbors + neighbors';
tee(fidLog, display_time_delay)

%% define all edges (taking care of neighborhoods with more than 2 areas)
%- based on the number of bridging watershed-line voxel between them (i.e. "neighbors" matrix)
tee(fidLog, 'Define edges... ')
tic
E = [];
ne2 = zeros(2,3);
for i=1:length(UN)
    if length(UN{i})>2
        UT = UN{i};
        nUT = length(UT);
        nTuple = nUT*(nUT-1)/2;
        ne = zeros(nTuple,3);
        c = 0;
        for j=1:nUT-1
            for k=j+1:nUT
                c=c+1;
                ne(c,:) = [UT(j), UT(k), neighbors(UT(j),UT(k))];
            end
        end
        %--- throw away tuples, which are not anymore touching after the merging operation
        if any(ne(:,3)==0)
            tee(fidLog, '\nWarning: removing lost borders from neighborhood\n')
            disp(ne)
            ne(ne(:,3)==0,:) = [];
            %disp(ne)
        end
        %--- disambiguate, only if more then two edges
        if size(ne,1)<=2
            UE = ne(:,1:2)';
        else
            
            [~,idx] = max(ne(:,3));
            UE = ne(idx,1:2);
            for j=1:nUT-2
                %--- test the left end
                idl = any(UE(1)==ne(:,1:2), 2);
                for k=2:length(UE)
                    idl = idl & ~any(UE(k)==ne(:,1:2), 2);
                end
                idx1 = find(idl);
                [~,idx2] = max(ne(idx1,3));
                ne2(1,:) = ne(idx1(idx2),:);
                %             UE(2+j) = temp(temp~=UE(end));
                %--- test the right end
                idl = any(UE(end)==ne(:,1:2), 2);
                for k=1:length(UE)-1
                    idl = idl & ~any(UE(k)==ne(:,1:2), 2);
                end
                idx1 = find(idl);
                [~,idx2] = max(ne(idx1,3));
                ne2(2,:) = ne(idx1(idx2),:);
                %             UE(2+j) = temp(temp~=UE(end));
                %---
                if ne2(1,3) > ne2(2,3)
                    UE = [    ne2(1, find(~ismember(ne2(1,1:2),UE))), UE];
                else
                    UE = [UE, ne2(2, find(~ismember(ne2(2,1:2),UE)))    ];
                end
            end
            UE = [UE(1:end-1); UE(2:end)];
            % UT(UE)
        end
        E = [E, UE];
    else
        E = [E, UN{i}];
    end
end
%---
E = unique(sort(E)','rows');
tee(fidLog, display_time_delay)



%% Remove end-nodes, which are mostly internal to another node
%- In particularly, check end-nodes, which are only one edge away from the next branching point
%- and remove the node if its internal surface is bigger than its external surface
tee(fidLog, 'Find end-nodes, which are mostly internal to other nodes... ')
tic
%--- find for each area the internal and external surface
Mtemp = single(ML) + (single(MS)*(MLn+1));
[internal, background] = find_clusters_neighbors( Mtemp , hdr.pixdim(2:4));
background(end) =[]; %--- remove the watershed-line area, which had the largest label
internal=internal(1:end-1, end);
%--- create graph and find end nodes of the graph
G = graph(E(:,1), E(:,2),[]);
deg = degree(G);
idxEnd = find(deg==1);
%--- find for each end node the connected node and if this has a degree >2 then test its surface exposure
labelOldNew = zeros(length(idxEnd),2);
idxEdgeRemove = zeros(length(idxEnd),1);
for i=1:length(idxEnd)
    idxEdge = find(any(E==idxEnd(i),2));
    edgeT = E(idxEdge,:);
    idxNext = edgeT(edgeT~=idxEnd(i));
    if deg(idxNext)>2 %--- test surface exposure
        if background(idxEnd(i))/internal(idxEnd(i))<=1
            labelOldNew(i, :) = [idxEnd(i), idxNext];
            idxEdgeRemove(i) = idxEdge;
        end
    end
end
labelOldNew(labelOldNew(:,1)==0, :) = [];
idxEdgeRemove(idxEdgeRemove==0) = [];

tee(fidLog, display_time_delay)
%% Remove such nodes and corresponding edges, if they exist

if ~isempty(idxEdgeRemove)
    
    %--- remove edges
    E(idxEdgeRemove,:) = [];
    
    
    %--- merge by changing the area labels
    tee(fidLog, 'Merge these nodes (%d) by changing the area labels... ', length(idxEdgeRemove))
    tic
    
    %--- change lable in variable U
    nChange = size(labelOldNew,1);
    for i=1:length(U)
        UT = U{i};
        for j=1:nChange
            UT(UT==labelOldNew(j,1)) = labelOldNew(j,2);
        end
        U{i} = unique(UT);
    end
    
    %--- change lable in variable UN
    nChange = size(labelOldNew,1);
    for i=1:length(UN)
        UT = UN{i};
        for j=1:nChange
            UT(UT==labelOldNew(j,1)) = labelOldNew(j,2);
        end
        UN{i} = unique(UT);
    end
    
    %--- get length of neighborhood tuples in U and UN
    Ul = cellfun(@length,U);
    UNl = cellfun(@length,UN);
    
    %--- give watershed-line voxel between merged areas an area-label (i.e. watershed-line voxel who are left with a single neighboring area)
    ML(Uind(Ul==1)) = [U{Ul==1}];
    
    %--- remove these watershed-line voxel with a single neighboring area from variables MS, U, Uind, Ul and in UN, UNl
    MS(Uind(Ul==1)) = 0;
    %---
    U(Ul==1,:) = [];
    Uind(Ul==1,:) = [];
    Ul(Ul==1) = [];
    %---
    UN(UNl==1,:) = [];
    UNl(UNl==1) = [];
    
    
    %--- change lable in variable ML
    for j=1:nChange
        ML(ML==labelOldNew(j,1)) = labelOldNew(j,2);
    end
    tee(fidLog, display_time_delay)
    
end

%% remove gaps in the list of labels, due to merging of areas
labelOld = unique(ML(:));
labelOld(labelOld==0) = []; %--- remove the zero, if present (should always be present as the background label)
%--- convert UN into a vector, for easier manipulation
UNl = cellfun(@length,UN);
UNlx = [0;UNl];
UNV = vertcat(UN{:});
%---
gapsN = max(labelOld) - length(labelOld);
if gapsN>0
    tee(fidLog, 'Making list of labels continuous, by removing %d missing labels... ', gapsN)
    tic
    labelOld = [0; labelOld(:)]; %--- add the "Zero", to allow detection of gaps of missing labels, starting with label "One"
    dd = diff(labelOld);
    idx = find(dd>1);
    offset = dd(idx) - 1;
    for i=length(idx):-1:1
        labelBorder = labelOld(idx(i));
        ML(ML>labelBorder) = ML(ML>labelBorder) - offset(i);
        E(E>labelBorder) = E(E>labelBorder) - offset(i);
        UNV(UNV>labelBorder) = UNV(UNV>labelBorder) - offset(i);
    end
    MLn = double(max(ML(:)));
    tee(fidLog, display_time_delay)
end
%--- convert the manipulated vector UNV back into the cell UN
for i=1:length(UN)
    UN{i} = UNV(sum(UNlx(1:i))+1:sum(UNlx(1:i+1)));
end


%% END of NODE and EDGE definitions


%% extract properties for each node
tee(fidLog, 'Extract node properties... ')
tic
T = regionprops('table',ML,{'Area','Centroid'});
T.cogI = T.Centroid(:,[2 1 3]); %--- x/y dimension are exchanged by "regionprops"
T.cog = (T.cogI .* repmat(hdr.pixdim(2:4), height(T), 1)) - repmat(hdr.pixdim(2:4)./2, height(T), 1);
%---
T.vol = T.Area .* prod(hdr.pixdim(2:4));
%---
tee(fidLog, display_time_delay)

%% calculate the distance between nodes, being connected by edges
tee(fidLog, 'Calculate edge length... ')
tic
D = zeros(size(E,1),1);
for i=1:size(E,1)
    D(i) = sum((T.cog(E(i,1),:) - T.cog(E(i,2),:)).^2)^.5;
end
tee(fidLog, display_time_delay)



%% CREATE THE GRAPH

%% create graph
G = graph(E(:,1), E(:,2),[], MLn);

%% degree of branching
T.degree = degree(G);

%% SAVE THE DATA
tee(fidLog, 'Saving data... \n')
tic
%--- determine the smallest data type allowing to save the volumetric data
typeStr = {'uint8','uint16','uint32','double'};
typeNr  = [2 4 8 62];
typeNeededStr = [];
typeNeededNr = [];
for i=1:3
    if MLn+1 <= intmax(typeStr{i})
        typeNeededStr = typeStr{i};
        typeNeededNr = typeNr(i);
        break
    end
end
if isempty(typeNeededStr)
    typeNeededStr = typeStr{4};
    typeNeededNr = typeNr(4);
end
tee(fidLog, 'Maximum area-label %d requires data type %s\n', MLn+1, typeNeededStr)

%--- save volume with area label + labeled watershed lines
MLWL = double(ML) + (double(MS)*(MLn+1));

%- create color map sorted by area size
MLP = regionprops('table', ML,'Area');
MLP.ColorIndex = MLP.Area;
MLP.ColorIndex = round(sqrt(MLP.ColorIndex)); %--- square root transformation
MLP.ColorIndex(MLP.ColorIndex>64) = 64; %--- 64^2=4096 corresponding to the upper limit (#voxel) where area volumes are clipped
cmapT = jet(64);
cmapT = cmapT(MLP.ColorIndex,:);
hdr.cmap = [0.5 0.5 0.5; cmapT; 1 1 1];

% fnameOut = [trunk,'_skeleton_watershed_test_end_nodes.nii'];
% fnameOut = [trunk,'_skeleton_watershed_areas.mat'];
fnameOut = [trunk,'_stack4_watershed_areas.mat'];
hdr = change_hdr_datatype_info(hdr,typeNeededNr);
write_3d_image_matrix(hdr,cast(MLWL, typeNeededStr),fnameOut);

%--- save the other data into a .mat file
% fnameOut = [trunk,'_skeleton_watershed_test_end_nodes.mat'];
% fnameOut = [trunk '_skeleton_watershed_graph.mat'];
fnameOut = [trunk '_stack4_watershed_graph.mat'];
save(fnameOut, 'T', 'E', 'D')
tee(fidLog, display_time_delay)

%%
tee(fidLog, 'Skeletonization completed!\n')
tee(fidLog, display_time_delay(hTimer))
fclose(fidLog);
fprintf('\n\n')

