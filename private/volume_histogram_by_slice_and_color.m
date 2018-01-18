function [h, idxSliceOk, prctD, histD, binEdges] = volume_histogram_by_slice_and_color(D, thrEffectiveResolution, mask, rangeBins, logTrans, hideFigure)
% function [h, idxSliceOk, prctD, histD, binEdges] = volume_histogram_by_slice_and_color(D, thrEffectiveResolution, mask, rangeBins, logTrans)
%
% This function plots the histogram for all slices in a 3D volume and for
% all color layers (max. 3 layers).
% In addition, it calculates the effective resolution of the image as the
% number of discrete levels used to encode the 1st-99th percentile range.
% It also returns the indices "idxSliceOk" of the first and last slice in the 3D volume
% having an effective resolution higher than "thrEffectiveResolution".
%% convert to logarithm?
if ~exist('logTrans','var') || isempty(logTrans)
    logTrans = true;
end

%% ignore threshold
if ~exist('thrEffectiveResolution','var') || isempty(thrEffectiveResolution)
    thrEffectiveResolution = -1;
end

%% does mask exist?
if exist('mask','var') && ~isempty(mask)
    maskExists = true;
    mask = mask>0;
    if ndims(D)==4 && ndims(mask)==3 %--- if mask is only 3D, assume that the mask is the same for all color layers
        mask = repmat(mask,1,1,1,size(D,4));
    elseif ndims(D)==4 && ndims(mask)==4 && size(D,4)>size(mask,4) %--- if mask is 4 D but has less color layers than D, assume that all missing color layers should be unmasked
        mask(:,:,:,size(mask,4)+1:size(D,4)) = true;
    end
else
    maskExists = false;
end

%% divide range of values into 256 bins (257 edges)
if exist('rangeBins','var')
    if isempty(rangeBins) %--- take zero to max range
        rangeBins = [0 max(D(:))];
    end
    interval = diff(rangeBins)/256;
    binEdges = rangeBins(1):interval:rangeBins(2);
    factVis = 256/diff(rangeBins);
else
    %--- assume that the range can be either [0,1] real numbers or [0,255] integers.
    if max(D(:))>1
        binEdges = 0:1:256;
        factVis = 1;
    else
        binEdges = 0:1/256:1;
        factVis = 255;
    end
end

%% calculate histograms and percentiles per z-slice and color-layer, for volume D
for iL = 1:size(D,4)
    for iS = 1:size(D,3)
        t = D(:,:,iS,iL);
        if maskExists
            t = t(mask(:,:,iS, iL));
        end
        histD(:,iS,iL) = histcounts(t(:),binEdges);
        % prctD(:,iS,iL) = prctile(t(:),0:5:100);
        prctD(:,iS,iL) = prctile(t(:),[1 99]);
    end
end
num = size(D,1)*size(D,2);
histD = histD./num*100; %--- convert in percentage

%% 
histMax = max(histD(:));

%% effective resolution -> number of intensity intervals used to encode the 1st-99th percentile range
effectiveResolution = squeeze(diff(prctD,[],1));

%% find longest contiguous range of slices with effective resolution above threshold
%--- find longest range for each color layer separately
for iL = 1:size(D,4)
    [idxSliceOkLayer(iL,1), idxSliceOkLayer(iL,2)] = find_longest_vector_segment(effectiveResolution(:,iL) > thrEffectiveResolution);
end
%--- find common range across color layer
idxSliceOk(1) = max(idxSliceOkLayer(:,1));
idxSliceOk(2) = min(idxSliceOkLayer(:,2));
%--- convert to double (needed for "text" function)
idxSliceOk = double(idxSliceOk);
%% plot histograms (as image with logarithmic color scale) and percentiles
%--- prepare figure
h(1) = figure;
if exist('hideFigure','var') && hideFigure
    h(1).Visible = 'off';
end
nL = size(D,4);
ScreenSz = get(0,'ScreenSize');
FigWidth = 556*nL;
h(1).Position = [ScreenSz(3)-FigWidth-7 591 FigWidth 754];
%---
colorLayer = {'Blue','Green','Red'};
for iL = 1:size(D,4)
    %surf(histDS(:,15:30,2))
    ha = subplot(1,nL,iL); hold on
    ha.Position = [(iL-1)*(1/nL)+0.06 0.07 (1/nL)-0.11 0.87];
    ha.Color = [0.7,0.7,0.7];
    ha.FontSize = 14;
    % title(sprintf('%s channel', colorLayer{iL}))
    % title(sprintf('%s channel: hist, prctil, hist-max', colorLayer{iL}))
    xlabel('slice in stack')
    ylabel('intensity')
    if logTrans
        title(sprintf('%s channel: log10(percentage(#pixel))', colorLayer{iL}))
        hi = imagesc(log10(histD(:,:,iL)),[-3 1.5]);
    else
        title(sprintf('%s channel: percentage(#pixel)', colorLayer{iL}))
        hi = imagesc(histD(:,:,iL),[0 ceil(histMax)]);
    end
    hi.AlphaData = (histD(:,:,iL)~=0)*1;
    % imagesc(histDS(:,:,iL),[0 15])
    axis tight
    colormap(jet(256))
    %--- add lines inidcating percentiles (every 10th percentile)
    % plot(prctD(1:2:end,:,iL)','-k','LineWidth',2)
    % plot(prctD(1:2:end,:,iL)',':w','LineWidth',2)
    plot(prctD(:,:,iL)' * factVis,'-r','LineWidth',2)
    plot(prctD(:,:,iL)' * factVis,':k','LineWidth',2)
    %--- add line indicating the maximum of the histograms at every slice
    %[~,histDmax]=max(histD(:,:,iL),[],1);
    %plot(histDmax,'-g','LineWidth',2);
    %plot(histDmax,':k','LineWidth',2);
    %--- add difference between 95th and 50th percentile
    % plot((prctDS(20,:,iL)-prctDS(11,:,iL)),'c','LineWidth',2)
    %--- add vertical lines indicating the longest range of good slices
    plot([idxSliceOkLayer(iL,1), idxSliceOkLayer(iL,1)], ylim, '--w', 'LineWidth', 2)
    plot([idxSliceOkLayer(iL,2), idxSliceOkLayer(iL,2)], ylim, '--w', 'LineWidth', 2)
    %---
    colorbar
    %--- adjust y-sclae
    ha.YTick = 1:255/5:256;
    ha.YTickLabel = round((ha.YTick-1)/factVis*100)/100;
end


%% plot effective resolution
if thrEffectiveResolution >= 0
    h(2) = figure;
    if exist('hideFigure','var') && hideFigure
        h(2).Visible = 'off';
    end
    ha = axes();
    ha.ColorOrder = [0 0 1; 0 1 0]; %--- blue, green
    ylim(binEdges([1,end]))
    hold on
    plot(effectiveResolution)
    title('effektive resolution')
    xlabel('slice')
    ylabel('percentile range [99th-1st]')
    ha.FontSize = 14;
    %--- find range of slices with good resolution
    hold on
    col='bgr';
    if ~all(isnan(idxSliceOk))
        for iL = 1:size(D,4)
            if any(isnan(idxSliceOkLayer(iL,:)))
                text(5, diff(ylim)*(1/(3-1*0.2)), 'WARNING: ALL SLICES BELOW THRESHOLD','FontSize',16, 'Color','r')
                text(5, diff(ylim)*(1/(3+3*0.2)), sprintf('         IN THE %s LAYER',colorLayer{iL}),'FontSize',16, 'Color',col(iL))
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
end
