function make_video_single_cell(h, trunk, ML, hdr, T, offset, MOther)

%%
% h = gcf;

%% set figure position and thereby the size
set(h,'Units', 'pixels')
% set(h,'Position', [600 200 1920 1080])
% set(h,'Position', [600 200 1496 1080])
set(h,'Position', [600 200 round(1496*(2/3)) round(1080*(2/3))])

%% create Video Writer Object
fname = [trunk, '_video.mp4'];
v = VideoWriter(fname, 'MPEG-4');
v.FrameRate = 20;
open(v)

%% find patches
hp = findall(h,'type','patch');
%--- change transparency
for j=1:length(hp)
    hp(j).FaceAlpha = 1;
end

%% find lights
lh = findall(h,'type','light');
camlight(lh,0,0)

%% take a first frame
drawnow nocallbacks
writeVideo(v,getframe(h));
% frame = getframe(h);
% writeVideo(v,frame);




%% sort areas by volume (i.e. cluster size) in descending order
% T = sortrows(T, 'Area', 'descend');
% %--- assign new label representing the new order
% % MLS = zeros(size(ML));
% MLS = ML;
% nL = height(T);
% for iL=1:nL
%     MLS(ML==T.Label(iL)) = iL;
% end
% %--- correct label for watershed line
% MLS(MLS==max(MLS(:))) = nL + 1;
% %--- add the other cells
% MLS(MOther>0) = nL + 2;
% %--- cluster sizes
% % nVox = T.Area;

%% reassign each area the volume of the area, to allow color coding according to area size
% T = sortrows(T, 'Area', 'descend');
T.ColorIndex = T.Area;
T.ColorIndex = round(sqrt(T.ColorIndex)); %--- square root transformation
T.ColorIndex(T.ColorIndex>64) = 64; %--- 64^2=4096 corresponding to the upper limit (#voxel) where ares volumes are clipped
%--- assign new label representing the new order
% MLS = zeros(size(ML));
MLS = ML;
nL = height(T);
for iL=1:nL
    MLS(ML==T.Label(iL)) = T.ColorIndex(iL);
end
%--- correct label for watershed line
MLS(ML==max(ML(:))) = 64 + 1;
%--- add the other cells
MLS(MOther>0) = 64 + 2;
%--- cluster sizes
% nVox = T.Area;
%% make sure no voxel with labels exeeding the color code are left
MLS(MLS>66) = 64 + 2;

%% create color map
cmap = [0.7 0.7 0.7; jet(64); 1 1 1; 0.3 0.3 0.3];

%% rotate the axes
tic
[az,el] = view;
% azp=0;
% el=14;
% el=0;
for i = 1:300
    %--- manipulate the view point
    %azp = azp + sin((-90+i)/180*pi)*1.5+2.5;
    %el = sin((i)/180*pi)*40;
    % azp = 179+i;
    % azp = -1+i;
    azp = i + az;
    % azp = -1+i/2 + az;
    %azp = -1+i;
    view(azp,el) %--- azimut and elevation
    camlight(lh,0,0)

    %--- get frame and save
    drawnow nocallbacks
    writeVideo(v,getframe(h));
end
toc

%% rotate the axes and increase transparency
tic
[az,el] = view;
% azp=0;
% el=14;
% el=0;
alpha = hp(1).FaceAlpha;
for i = 1:60+360+260
    %--- manipulate the view point
    %azp = azp + sin((-90+i)/180*pi)*1.5+2.5;
    %el = sin((i)/180*pi)*40;
    % azp = 179+i;
    % azp = -1+i;
    azp = i + az;
    % azp = -1+i/2 + az;
    %azp = -1+i;
    view(azp,el) %--- azimut and elevation
    camlight(lh,0,0)
    %--- change transparency
    if alpha>0.4 && mod(i-1,6)==0
        alpha = alpha - 0.01;
        for j=1:length(hp)
            hp(j).FaceAlpha = alpha;
        end
    end
    %--- get frame and save
    drawnow nocallbacks
    writeVideo(v,getframe(h));
end
toc

%% rotate the axes and increase transparency
tic
[az,el] = view;
% azp=0;
% el=14;
% el=0;
for i = 1:100
    %--- manipulate the view point
    %azp = azp + sin((-90+i)/180*pi)*1.5+2.5;
    %el = sin((i)/180*pi)*40;
    % azp = 179+i;
    % azp = -1+i;
    azp = i + az;
    % azp = -1+i/2 + az;
    %azp = -1+i;
    view(azp,el) %--- azimut and elevation
    camlight(lh,0,0)
    %--- change transparency
    if alpha>0.1 && mod(i-1,6)==0
        alpha = alpha - 0.01;
        for j=1:length(hp)
            hp(j).FaceAlpha = alpha;
        end
    end
    %--- get frame and save
    drawnow nocallbacks
    writeVideo(v,getframe(h));
end
toc

%% make surfaces totaly transparent
% for j=1:length(hp)
%     hp(j).FaceAlpha = 0;
% end
%% take a few still frames before showing the slices
drawnow nocallbacks
frame = getframe(h);
for i=1:20
    writeVideo(v,frame);
end

%--- change transparency
for j=1:length(hp)
    hp(j).FaceAlpha = 0.2;
end

%% create slices
%--- delete camplight (makes surface appear flat)
delete(lh)
%---
% Mtemp = single(MLS) + (single(MS)*(nMLS+1));
% Ms = single(Mtemp);
MLS = flip(rot90(MLS),1);
XV = ((1:size(MLS,1)) +  offset(2)) * hdr.pixdim(2) - hdr.pixdim(2)/2;
YV = ((1:size(MLS,2)) +  offset(1)) * hdr.pixdim(3) - hdr.pixdim(3)/2;
ZV = ((1:size(MLS,3)) +  offset(3)) * hdr.pixdim(4) - hdr.pixdim(4)/2;
[X,Y,Z] = meshgrid(YV,XV,ZV);
%---

colormap(cmap)

%--- walk through the volume in z-direction
xs = [];
ys = [];
zs = [];
hs = [];
zlimOrig  = zlim();
for i=size(MLS,3):-1:5
    zs = ZV(i-4);
    delete(hs);
    hs = slice(X,Y,Z, MLS,xs,ys,zs);
    for j=1:length(hs)
        hs(j).EdgeAlpha =0;
        hs(j).FaceAlpha =1;
    end
    if i<size(MLS,3)
        zlim([zlimOrig(1) ZV(i)])
    end
    drawnow nocallbacks
    frame = getframe(h);
    for j=1:7
        writeVideo(v,frame);
    end
%     pause(0.5)
end

%% rotate the axes and decrease transparency
% tic
% [az,el] = view;
% % azp=0;
% % el=14;
% % el=0;
% alpha = hp(1).FaceAlpha;
% for i = 1:360
%     %--- manipulate the view point
%     %azp = azp + sin((-90+i)/180*pi)*1.5+2.5;
%     %el = sin((i)/180*pi)*40;
%     % azp = 179+i;
%     % azp = -1+i;
%     azp = i + az;
%     % azp = -1+i/2 + az;
%     %azp = -1+i;
%     view(azp,el) %--- azimut and elevation
%     camlight(lh,0,0)
%     %--- change transparency
%     if i<= 0.3*360 && mod(i,6)==0
%         for j=1:length(hp)
%             hp(j).FaceAlpha = alpha + i/360;
%         end
%     end
%     drawnow nocallbacks
%     %--- get frame and save 
%     writeVideo(v,getframe(h));
% end
% toc







%% pan the axes along the x-axis
% ha = findall(h,'type','axes');
% %---
% xlimFull = xlim(ha);
% %---
% tic
% for i=xlimFull(1):0.1:xlimFull(2)-1
%     xlim(ha, [i xlimFull(2)])
%     %---
%     drawnow nocallbacks
%     % pause(0.5)
%     writeVideo(v,getframe(h));
% end
% toc


%% save movie to file

close(v)

