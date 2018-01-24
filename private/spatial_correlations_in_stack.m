function [h, idxSliceOk, spatialCorrSmoothed] = spatial_correlations_in_stack(img1, direction, thr, hideFigure)
% function [h, idxSliceOk] = spatial_correlations_in_stack(img1, direction, thr)
nL = size(img1,4);
nS = size(img1,3);
if strcmp(direction,'vertical')
    %--- Calculate spatial correlations between succesive slices
    spatialCorr = zeros(nS-1,nL);
    for iL =1:nL
        for iS = 1:nS-1
            %         r = corrcoef(img1(:,:,iS+1,iL), img1(:,:,iS,iL));
            %         spatialCorr(iS,iL) = r(1,2);
            spatialCorr(iS,iL) = corr2(img1(:,:,iS+1,iL), img1(:,:,iS,iL));
            
        end
    end
elseif strcmp(direction,'horizontal')
    %--- Calculate spatial correlations within slice
    spatialCorrWithin = zeros(4, nS,nL);
    for iL =1:nL
        for iS = 1:nS
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
    if nL==1
        spatialCorr = spatialCorr';
    end
end
%% smooth across slices
spatialCorrSmoothed = nan(size(spatialCorr));
for iL =1:nL
    spatialCorrSmoothed(:,iL) = smooth_moving_average(spatialCorr(:,iL), 2);
end
%% find longest contiguous range of slices with spatial correlations above threshold
%--- find longest range for each color layer separately
for iL = 1:nL
    [idxSliceOkLayer(iL,1), idxSliceOkLayer(iL,2)] = find_longest_vector_segment(spatialCorrSmoothed(:,iL) > thr);
end
%--- find common range across color layer
idxSliceOk(1) = max(idxSliceOkLayer(:,1));
idxSliceOk(2) = min(idxSliceOkLayer(:,2));
%--- convert to double (needed for "text" function)
idxSliceOk = double(idxSliceOk);
%%
col='bgr';
layerName = {'BLUE','GREEN','RED'};
%--- plot the correlations
h = figure;
if exist('hideFigure','var') && hideFigure
    h.Visible = 'off';
end
ha = axes(); hold on
for iL = 1:nL
    plot(spatialCorrSmoothed(:,iL),'k')
    plot(spatialCorr(:,iL), col(iL), 'LineWidth',2)
end
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
if ~all(isnan(idxSliceOk))
    for iL = 1:nL
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
