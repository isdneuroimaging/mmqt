function [ cl, voxN, voxType, volume, area ] = find_individual_clusters_find_surface_anisotrop( volT, thr, pixdim )
% function [ cl, voxN, voxType, volume, area ] = find_individual_clusters_find_surface_anisotrop( volT, thr, pixdim )
% written by Benno, 01.04.2014
%%%--- General notes:
%--- The algorithm loops systematically through the volume of voxels and
%--- tests each voxel against the threshold. As soon as a suprathreshold voxel
%--- is found, it starts a see-growing algorithm. This algorithm tests each neighboring
%--- voxel against the threshold. If a neighboring voxel is supra-threshold,
%--- then it gets the same cluster ID as the seed-voxel and is added to a
%--- queue. The seed-growing algorithm walks through this queue and repeats for
%--- each element the check of its neighbors, until the end of this queue is
%--- reached. In addition to testing against each neighboring voxel against the
%--- threshold "thr", the seed-growing algorithm checks whether a neighbor is
%--- already part of a cluster (i.e. whether a cluster ID was already
%--- atributed to it). 
%--- Only voxels bordering with a side-wall to the seed-voxel are considered neighbors. 
%--- This means that clusters, which are bordering only over a corner or
%--- edge against each other are not merged together.

%% Initialize some variables
%--- Add an aditional layer of "zero"-voxels around the input volume "volT",
%--- in order to allow the seed growing algorithm to continue woring at
%--- the edges of the original input volume, without additional checks for
%--- the location of the voxels within the volume. These voxels have to be
%--- removed again, at the end.
% vol = zeros(size(volT)+2,class(volT));
vol = zeros(size(volT)+2,class(volT));
vol(2:end-1,2:end-1,2:end-1) = volT;

%--- Create an empty volume, to store information of the identity of found clusters
cl = uint16(zeros(size(vol)));
voxType = uint16(zeros(size(vol))); %--- only needed in case that you want to test whether a voxel is a surface voxel of a cluster
%--- Define additional variables for...
clN = 0; %--- number of found clusters
voxN = []; %--- number of voxels per cluster; this vector grows with each newly found cluster
volume = []; %--- volume per cluster
area = []; %--- area per cluster
pixelArea(1) = prod(pixdim([2,3]));
pixelArea(2) = prod(pixdim([1,3]));
pixelArea(3) = prod(pixdim([1,2]));
voxelVolume = prod(pixdim);
%--- Create a queue (2D matrix) for the seed growing algorithm;
%--- NOTE: The length of the queue limits how big a cluster can maximally grow.
%---       If the end of the queue is reached, then the algorithm is
%---       aborted.
q=int16(zeros(50000000,3));

%% go through brain-volume and grow a seed for each encounter of a superthreshold voxel not already being classified
% tic
for i =2:size(vol,1)-1
    for j =2:size(vol,2)-1
        for k =2:size(vol,3)-1
            if vol(i,j,k)>thr && cl(i,j,k)==0 %--- this means that a new cluster was encountered! Grow the full cluster from this seed! 
                %--- increment counter of found clusters
                clN = clN+1;
                hh = zeros(1,3);
                % area(clN) = 0;
                %--- add the voxel to the queue, in order to check its neighbors
                q(1,1)=i; q(1,2)=j; q(1,3)=k;
                %--- assign the cluster ID to the voxel
                cl(i,j,k)=clN;
                m=1; %--- initalize the counter for the while loop, which loops through the queue "q", until the most recently added element (voxel) was processed
                n=1; %--- index of the most recently added element (voxel) in the queue "q"; with other words: "n" counts how many voxels are currently in the growing cluster.
                while m<=n
                    x=q(m,1); y=q(m,2); z=q(m,3);
%                     if cl(x, y, z) %--- check weather already tested, in case that the voxel was inserted several times.
%                         m=m+1;
%                     else
                        % h=0; %--- only needed in case that you want to test whether a voxel is a surface voxel of a cluster
                        h = zeros(1,3);
                        %--- create a matrix with the coordinates of the neighboring voxels
                        cc=[x-1, y, z;...
                            x+1, y, z;...
                            x, y-1, z;...
                            x, y+1, z;...
                            x, y, z-1;...
                            x, y, z+1];
                        for g=1:6 %--- test all 6 neighboring voxels against the threshold
                            %--- only if not already tested AND if supra-threshold, then insert the voxel into the queue! This ensures, that no voxel will be added twice to the queue!
                            if vol(cc(g,1),cc(g,2),cc(g,3))>thr %--- check whether the voxel is suprathreshold
                                % h=h+1; %--- only needed in case that you want to test whether a voxel is a surface voxel of a cluster
                                if ~cl(cc(g,1),cc(g,2),cc(g,3)) %--- check whether voxel was already added to the cluster
                                    %--- assign the cluster ID to the voxel
                                    cl(cc(g,1),cc(g,2),cc(g,3))=clN;
                                    %--- add the voxel to the queue, in order to check its neighbors
                                    n=n+1;
                                    q(n,1)=cc(g,1); q(n,2)=cc(g,2); q(n,3)=cc(g,3);
                                end
                            else
                                idxh = round(g/2);
                                h(idxh) = h(idxh) + 1;
                            end
                        end
                        %--- check whether it is a surface or core voxel (if needed, then assignments to variables "voxType" and "h" have to be uncommented in previous rows)
                        % if h==6; voxType(x,y,z)=2; else voxType(x,y,z)=1; end
                        if sum(h)==0; voxType(x,y,z)=2; else voxType(x,y,z)=1; end %--- sum(h)>0 means that it is a surface voxel, exposed to the surface at least with one of it's walls
                        hh = hh + h; %--- sum of voxel-walls exposed to surface
                        % area(clN) = area(clN) + (6-h);
                        %--- increment the while-loop counter
                        m=m+1;
                        if n>size(q,1)-6; fprintf('WARNING: cluster aborted at %d voxels\n', n)
                            break
                        end
%                     end
                end
                voxN(clN) = n;
                area(clN) = sum(hh .* pixelArea);
                volume(clN) = n * voxelVolume;
                % fprintf('%d\n',n)
            end
        end
    end
end
% toc
%% go through brain-volume and renumber clusters according to their size
%--- remove extra layer of voxels
cl=cl(2:end-1,2:end-1,2:end-1);
voxType=voxType(2:end-1,2:end-1,2:end-1); %--- only needed in case that you want to test whether a voxel is a surface voxel of a cluster
%--- sort voxN and get the index
[voxN, idx]=sort(voxN, 'descend');
area = area(idx);
volume = volume(idx);
[~, idx]=sort(idx);
% tic
for x =1:size(cl,1)
    for y =1:size(cl,2)
        for z =1:size(cl,3)
            if cl(x,y,z)>0
                %--- exchange it's ID
                cl(x,y,z) = idx(cl(x,y,z));
            end
        end
    end
end
% toc
end

