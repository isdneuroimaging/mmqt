function adjust_viewing_options(ha, hdr, T)
%%
ha.FontSize = 14;
%% viewing options
view(3)
% view(-38,23)
% view(0,90)
axis vis3d;
axis tight
axis equal
rotate3d on

xlabel('X [\mum]')
ylabel('Y [\mum]')
zlabel('Z [\mum]')
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% dimRealWorld = sizeVolume .* hdr.pixdim(2:4);
% xlim([0, dimRealWorld(1)])
% ylim([0, dimRealWorld(2)])
% zlim([0, dimRealWorld(3)])
% axis off
% ha.Position = [0.2 0.2 0.6 0.6]
% ha.Position = [0.05 0.05 0.9 0.9];
%%
limits = [min(T.cog)-5;...
    max(T.cog)+5]';
try
    xlim([limits(1,1), limits(1,2)]);
    ylim([limits(2,1), limits(2,2)])
    zlim([limits(3,1), limits(3,2)])
catch
end
