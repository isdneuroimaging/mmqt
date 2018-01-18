function [neighbors, background] = find_clusters_neighbors( L, pixdim )
% function [neighbors, background] = find_clusters_neighbors( L, pixdim )
% written by Benno, 21.01.2017
%
% Input:
%       L - image matrix with labeled areas (2D or 3D are possible)
%       pixdim - (optional) voxel dimensions; vector of length 3, giving edge length in x, y, and z direction;
%                if given, then the actual area of contact between
%                neighboring areas will be calculated
%
% Output:
%       neighbors - symmetric matrix, containing for each pair of areas the number of voxel-faces over which these areas touch each other
%       background - vector containing the number of voxel-faces over which each area touches the background
%
%%%--- General notes:
%Script loops through all voxel in the picture and tests the identity of 
%the six all neighboring voxel touching with a surface.
%The result will be the matrix "neighbors" showing the number of surfaces over which
%two clusters are in contact with each other. The diagonal of the matrix
%will be the number of internal surfaces of the corresponding cluster
%multiplied by two (because each surface is tested twice)
%Voxel with label Zero are considered background. Touches with the
%background are recorded in the separate output "background"


%% Define the voxel-neighborhood to be tested for each voxel (accounting for the border of the image)
%--- define the coordinates of the neighboring voxel to be tested, relative to the center voxle coordinates
neCoord = [-1 0 0;...
             1 0 0;...
             0 -1 0;...
             0 1 0;...
             0 0 -1;...
             0 0 1];
%--- surface area of voxel faces
if exist('pixdim','var') && ~isempty(pixdim)
    area = [pixdim(2)*pixdim(3), pixdim(1)*pixdim(3), pixdim(1)*pixdim(2)];
    area = reshape([area; area],[],1);
end
%--- get indices and subscripts for all voxel in the volume
% ind = 1:numel(L);
%--- get indices and subscripts for all non-Zero voxel in the volume
ind = find(L>0);
[x,y,z] = ind2sub(size(L), ind);
%--- define for each voxel, which neighbors to test (i.e. ne2test)
ne2test = false(length(ind),6);
ne2test(:,1) =  x>1;
ne2test(:,2) =  x<size(L,1);
ne2test(:,3) =  y>1;
ne2test(:,4) =  y<size(L,2);
ne2test(:,5) =  z>1;
ne2test(:,6) =  z<size(L,3);


%% loop through all voxel in the picture and test identity of all touching voxel
%--- increase the counter of the matrix cell at the row and column
%--- corresponding to the identity of the respectivly tested seed and neighbor
nCl = max(L(:));
neighbors = zeros(nCl,nCl+1); %--- one additional column for the background (i.e. label==0)
if exist('pixdim','var') && ~isempty(pixdim) %--- in case that the pixeldimensions are known calculate the exact surface area
    for i =1:length(ind)
        label = L(ind(i));
        cc = repmat([x(i) y(i) z(i)],sum(ne2test(i,:)),1) + neCoord(ne2test(i,:),:);
        areaT = area(ne2test(i,:));
        for g=1:size(cc,1) %--- test all neighboring voxels against the threshold
            labelNeighbor = L(cc(g,1),cc(g,2),cc(g,3));
            neighbors(label, labelNeighbor+1) = neighbors(label, labelNeighbor+1) + areaT(g);
        end
    end
else
    for i =1:length(ind)
        label = L(ind(i));
        cc = repmat([x(i) y(i) z(i)],sum(ne2test(i,:)),1) + neCoord(ne2test(i,:),:);
        for g=1:size(cc,1) %--- test all neighboring voxels against the threshold
            labelNeighbor = L(cc(g,1),cc(g,2),cc(g,3));
            neighbors(label, labelNeighbor+1) = neighbors(label, labelNeighbor+1) + 1;
        end
    end
end
%% separate the column with the background
background = neighbors(:,1); %--- first column was the contacts with the background (i.e. label==0)
neighbors = neighbors(:,2:end);

%% OLD VERSION OF: 
% %% Initialize some variables
% %--- Add an aditional layer of "zero"-voxels around the input volume "volT",
% %--- in order to allow the seed growing algorithm to continue woring at
% %--- the edges of the original input volume, without additional checks for
% %--- the location of the voxels within the volume. These voxels have to be
% %--- removed again, at the end.
% sz = size(L);
% if length(sz)==2
%     sz(3) = 1;
% end
% L = zeros(sz+2,class(L));
% L(2:end-1,2:end-1,2:end-1) = L;
% %% loop through all voxel in the picture and test identity of all touching voxel
% %--- increase the counter of the matrix cell at the row and column
% %--- corresponding to the identity of the respectivly tested seed and neighbor
% % tic
% nCl = max(L(:));
% neighbors = zeros(nCl,nCl);
% background = zeros(nCl,1);
% for x =2:size(L,1)-1
%     for y =2:size(L,2)-1
%         for z =2:size(L,3)-1
%             idx = L(x,y,z);
%             if idx>0
%                 cc=[x-1, y, z;...
%                     x+1, y, z;...
%                     x, y-1, z;...
%                     x, y+1, z;...
%                     x, y, z-1;...
%                     x, y, z+1];
%                 for g=1:6 %--- test all 6 neighboring voxels against the threshold
%                     idxN = L(cc(g,1),cc(g,2),cc(g,3));
%                     if idxN == 0
%                         background(idx) = background(idx) + 1;
%                     else
%                         neighbors(idx, idxN) = neighbors(idx, idxN)+1;
%                     end
%                 end
%             end
%         end
%     end
% end
