function [h, idxSliceOk, spatialCorrSmoothed] = spatial_correlations_in_stack(img1, direction, thr, hideFigure)
% function [h, idxSliceOk] = spatial_correlations_in_stack(img1, direction, thr)

if strcmp(direction,'vertical')
    %--- Calculate spatial correlations between succesive slices
    spatialCorr = zeros(size(img1,3)-1,size(img1,4));
    for iL =1:2
        for iS = 1:size(img1, 3)-1
            %         r = corrcoef(img1(:,:,iS+1,iL), img1(:,:,iS,iL));
            %         spatialCorr(iS,iL) = r(1,2);
            spatialCorr(iS,iL) = corr2(img1(:,:,iS+1,iL), img1(:,:,iS,iL));
            
        end
    end
elseif strcmp(direction,'horizontal')
    %--- Calculate spatial correlations within slice
    spatialCorrWithin = zeros(4, size(img1,3),size(img1,4));
    for iL =1:2
        for iS = 1:size(img1, 3)
            %         r = corrcoef(img1(:,:,iS+1,iL), img1(:,:,iS,iL));
            %         spatialCorr(iS,iL) = r(1,2);
            spatialCorrWithin(1, iS,iL) = corr2(img1(1:end-1,1:end-1,iS,iL), img1(2:end,2:end,iS,iL)); %--- 1. diagonal
            spatialCorrWithin(2, iS,iL) = corr2(img1(2:end,1:end-1,iS,iL), img1(1:end-1,2:end,iS,iL)); %--- 2. diagonal
            spatialCorrWithin(3, iS,iL) = corr2(img1(1:end-1,:,iS,iL), img1(2:end,:,iS,iL)); %--- horizontal neighbors
            spatialCorrWithin(4, iS,iL) = corr2(img1(:,1:end-1,iS,iL), img1(:,2:end,iS,iL)); %--- vertical neighbors
        end
    end
    %--- average across directions
    spatialCorr = squeeze(mean(spatialCorrWithin));
end
%% smooth across slices
spatialCorrSmoothed = nan(size(spatialCorr));
for iL =1:size(spatialCorr,2)
    spatialCorrSmoothed(:,iL) = smooth_moving_average(spatialCorr(:,iL), 2);
end
%% find longest contiguous range of slices with spatial correlations above threshold
%--- find longest range for each color layer separately
for iL = 1:size(img1,4)
    [idxSliceOkLayer(iL,1), idxSliceOkLayer(iL,2)] = find_longest_vector_segment(spatialCorrSmoothed(:,iL) > thr);
end
%--- find common range across color layer
idxSliceOk(1) = max(idxSliceOkLayer(:,1));
idxSliceOk(2) = min(idxSliceOkLayer(:,2));
%--- convert to double (needed for "text" function)
idxSliceOk = double(idxSliceOk);
%%
%--- plot the correlations
h = figure;
if exist('hideFigure','var') && hideFigure
    h.Visible = 'off';
end
ha = axes();
plot(spatialCorrSmoothed(:,1),'k')
hold on
plot(spatialCorrSmoothed(:,2),'k')
plot(spatialCorr(:,1),'b','LineWidth',2)
plot(spatialCorr(:,2),'g','LineWidth',2)
ylim([0 1])
if strcmp(direction,'vertical')
    title('spatial correlations between succesive slices')
else
    title('spatial correlations within slice')
end
xlabel('slice')
ylabel('Pearson''s r')
% ylabel('corr 2D')
ha.FontSize = 14;
%--- find range of slices with good resolution
hold on
col='bgr';
layerName = {'BLUE','GREEN','RED'};
if ~all(isnan(idxSliceOk))
    for iL = 1:size(img1,4)
        if any(isnan(idxSliceOkLayer(iL,:)))
            text(5, diff(ylim)*(1/(3-1*0.2)), 'WARNING: ALL SLICES BELOW THRESHOLD','FontSize',16, 'Color','r')
            text(5, diff(ylim)*(1/(3+3*0.2)), sprintf('         IN THE %s LAYER',layerName{iL}),'FontSize',16, 'Color','r')
        else
            % [idxSliceOkLayer(iL,1), idxSliceOkLayer(iL,2)] = find_longest_vector_segment(effectiveResolution(:,iL) > thrEffectiveResolution);
            plot([idxSliceOkLayer(iL,1), idxSliceOkLayer(iL,1)], ylim, col(iL), 'LineWidth', 2)
            plot([idxSliceOkLayer(iL,2), idxSliceOkLayer(iL,2)], ylim, col(iL), 'LineWidth', 2)
        end
    end
end
%--- find common window with good signal-to-noise ratio in both layers
% idxSliceOk(1) = max(idxSliceOkLayer(:,1));
% idxSliceOk(2) = min(idxSliceOkLayer(:,2));
if all(isnan(idxSliceOk))
    text(5, diff(ylim)*(2/(3+1*0.2)), 'WARNING: ALL SLICES BELOW THRESHOLD','FontSize',16, 'Color','r')
else
    plot([idxSliceOk(1), idxSliceOk(1)], ylim, ':k', 'LineWidth', 2)
    plot([idxSliceOk(2), idxSliceOk(2)], ylim, ':k', 'LineWidth', 2)
    text(idxSliceOk(1), diff(ylim)*(2/(3+1*0.2)), num2str(idxSliceOk(1)),'FontSize',16)
    text(idxSliceOk(2), diff(ylim)*(2/(3+2*0.2)), num2str(idxSliceOk(2)),'FontSize',16)
end
