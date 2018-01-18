function [hl, hp, he] = add_isosurface(M, hdr, offset, cellTouches, addLight, alpha)
%% add isosurface

%%
if ~exist('alpha','var') || isempty(alpha)
    alpha = 0.4;
end
%%
whichplane = {'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'};
cellTouches = reshape(cellTouches',6,1);
%%
%--- change x and y dimension 
pixdim = hdr.pixdim([3 2 4]);
offset = offset([2 1 3]);
M90 = flip(rot90(M,1),1);
%---
% sigma = 0.22 ./pixdim;
sigma = 0.15 ./pixdim;
M90s = imgaussfilt3(single(M90),sigma);
% M90s = M90;
%---
% figure
% [X,Y,Z] = meshgrid((1:size(M90,1))*pixdim(1), (1:size(M90,2))*pixdim(2), (1:size(M90,3))*pixdim(3));
vx = ((1:size(M90,2)) + offset(2)) *pixdim(2);
vy = ((1:size(M90,1)) + offset(1)) *pixdim(1);
vz = ((1:size(M90,3)) + offset(3)) *pixdim(3);
[X,Y,Z] = meshgrid(vx, vy, vz);
% patch(isocaps(X,Y,Z,M90s,0.01), 'FaceColor', [0.5 0 1], 'EdgeColor', 'none')
for i=1:length(whichplane)
    if cellTouches(i)
        he(i) = patch(isocaps(X,Y,Z,M90s,0.1, whichplane{i}), 'FaceColor', [1 0 0.5], 'EdgeColor', 'none', 'BackFaceLighting', 'unlit');
    else
        he(i) = patch(isocaps(X,Y,Z,M90s,0.1, whichplane{i}), 'FaceColor', [0.6 0.6 1], 'EdgeColor', 'none', 'BackFaceLighting', 'unlit');
    end
end
hp = patch(isosurface(X,Y,Z,M90s,0.1));
%---
% hp.FaceColor = [0.5 0.5 1];
hp.FaceColor = [0.6 0.6 1];
hp.EdgeColor = 'none';
hp.FaceAlpha = alpha;
hp.SpecularStrength = 0.4;
hp.BackFaceLighting = 'unlit';
%---
if ~exist('addLight', 'var') || addLight
    hl = camlight(-25, 25);
    lighting gouraud
else
    hl = [];
end
%---
ylim([vy(1)-1 vy(end)+1])
xlim([vx(1)-1 vx(end)+1])
zlim([vz(1)-1 vz(end)+1])
