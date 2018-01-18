function mmqt_segment_image(trunk, figureScaling, figureHide)
% function mmqt_segment_image(fnameCZI, figureScaling, figureHide)
% function mmqt_segment_image(fnameMat, figureScaling, figureHide)
%
% Segment foreground from background. In particular, cell nuclei and
% micgroglia will be segmented from the coresponding color layer.
% Input should be a Z-stack with two color layers:
%   1. Layer: DAPI staining of nuclei
%   2. Layer: anti-Iba1 steining of microglia
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
fprintf('\n\nSegmenting image:\n\n')

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

%% Specify some parameters
rSoma = 2.2;
thrVolume = 10;

%% Create log file
hTimer = tic;
fnameLog = [trunk, '_log1_segmentation.txt'];
fid = fopen(fnameLog, 'w');
fprintf(fid, '%s\n', date);
fclose(fid);

%% Reconstruct MAT file name with raw data
fname = [trunk '_stack1_raw.mat'];
if ~exist(fname, 'file')
    fname = [trunk '.mat'];
    if ~exist(fname, 'file')
        error('%s\n%s\n  %s\n  %s\n', 'File with raw data Z-stack not found!', 'The file-name was expected to be either one of these:', [trunk '_stack1_raw.mat'], [trunk '.mat'])
    end
end
%% Read MAT file
[hdr, img1] = read_3D_image_matrix(fname);
load(fname, 'img', 'hdr');

%% normalize to range [0,1]
img2 = img1./max(img1(:));


%% smooth the data
tee(fnameLog, 'Smoothing the data slice by slice... \n')
imgS = zeros(size(img2));
sigma = 0.3 ./ hdr.pixdim(2:3); %--- pixel dimension is given in micrometer
%--- display in command window
tee(fnameLog, '- Voxel dimension: %.3f %.3f %.3f\n', hdr.pixdim(2:4))
tee(fnameLog, '- Sigma:           %g %g\n', sigma)
%--- smooth the data
tic
for iL=1:size(img2,4)
    %--- slice by slice
    for iS=1:size(img2,3)
        imgS(:,:,iS,iL) = imgaussfilt(img2(:,:,iS,iL),sigma);
    end
end
tee(fnameLog, display_time_delay)


%% check out histograms and effective resolution (number of discrete levels used to encode the 1st-99th percentile range)
[hf, ~, prctD, histD] = volume_histogram_by_slice_and_color(img1, 45, [], [], [], figureHide);

%% save picture of histogram
%--- histogram
F = getframe(hf(1));
imwrite(F.cdata, [trunk,'_prepro2_histograms_raw.png'])
%--- effective resolution
F = getframe(hf(2));
imwrite(F.cdata, [trunk,'_prepro1_effective_resolution.png'])

%% spatial correlations
thrSpatialCorr = 0.78;
tee(fnameLog, 'Calculate spatial correlations ...  ')
tic
%--- Z direction (correlations between succesive slices)
%- Note: spatialCorrSmoothed is not used, but saved in  a mat-file
[h, idxSliceOk, spatialCorrSmoothed] = spatial_correlations_in_stack(img2, 'vertical', thrSpatialCorr, figureHide);
%--- save picture of spatial correlations (across slice)
F = getframe(h);
imwrite(F.cdata, [trunk,'_prepro1_spatial_corrZ.png'])

%--- XY direction (within slice)
[h] = spatial_correlations_in_stack(img2, 'horizontal', thrSpatialCorr, figureHide);
%--- save picture of spatial correlations (within slice)
F = getframe(h);
imwrite(F.cdata, [trunk,'_prepro1_spatial_corrXY.png'])
tee(fnameLog, display_time_delay)
%%
if any(isnan(idxSliceOk))
    warning('None of the slices has the required quality, according to spatial correlations!')
    warning('Analysing the whole stack!')
    idxSliceOk = [1 size(img1,3)];
    %--- write warning to log file
    fid = fopen([trunk,'_warnings.txt'], 'a');
    fprintf(fid, 'None of the slices has the required quality, according to spatial correlations!\n');
    fprintf(fid, ' => Analysing the whole stack!\n');
end

%% save result of spatial correlations in z-direction
save([trunk,'_prepro1_spatial_corrZ_sliceOk.mat'], 'idxSliceOk', 'spatialCorrSmoothed', 'thrSpatialCorr')

%% equalize histograms
%--- define a norm histogram to which all histograms shell be matched using histeq.m
for iL=1:size(prctD,3)
    %--- as norm, select percentiles from slice with largest difference between the 99th and 1th percentile
    prctRange = prctD(2,:,iL)-prctD(1,:,iL);
    prctRmax = max(prctRange);
    %--- consider also slices where the criterion has a value of at least 95% of the maximum
    idx = find(prctRange>prctRmax*0.95);
    %--- calculate the norm histogram as the mean histogram across the selected slices
    histNorm(:,iL) = mean(histD(:,idx,iL),2);
end

%--- equalize histograms using "histeq"
imgEq = zeros(size(img2));
for iL = 1:size(img2,4)
    for iS = 1:size(img2,3)
        imgEq(:,:,iS,iL) = histeq(img2(:,:,iS,iL),histNorm(:,iL));
    end
end

%% plot histogram
hf = volume_histogram_by_slice_and_color(imgEq, [], [], [], [], figureHide);
%%
F = getframe(hf);
imwrite(F.cdata, [trunk,'_prepro3_histograms_equalized.png'])

%% normalize each color layer to a certain percentile range (for visualization)
prct = prctile(reshape(imgEq(:,:,:,1),[],1),[1 99]);
imgEq(:,:,:,1) = imgEq(:,:,:,1) - prct(1);
imgEq(:,:,:,1) = imgEq(:,:,:,1)./diff(prct);
prct = prctile(reshape(imgEq(:,:,:,2),[],1),[1 99]);
imgEq(:,:,:,2) = imgEq(:,:,:,2) - prct(1);
imgEq(:,:,:,2) = imgEq(:,:,:,2)./diff(prct);
imgEq(imgEq<0) = 0;
imgEq(imgEq>1) = 1;


%% show the whole stack
h = show_cells_3D(imgEq,hdr,[],[],figureScaling, [], figureHide);
%--- save picture
F = getframe(h);
imwrite(F.cdata, [trunk,'_ortho2_equalized.png']);
%% crop the stack
imgEqCrop = imgEq(:,:,idxSliceOk(1):idxSliceOk(2),:);
imgSCrop = imgS(:,:,idxSliceOk(1):idxSliceOk(2),:);
% imgSCrop = imgS;

%% show the cropped stack
h = show_cells_3D(imgEqCrop,hdr,[],[],figureScaling, [], figureHide);
%--- save picture
F = getframe(h);
imwrite(F.cdata, [trunk,'_ortho3_equalizedCrop.png']);

%% Determine threshold for each slice separately
tee(fnameLog, 'Determine thresholds for each slice... \n')
tic
hf = figure;
if figureHide
    hf.Visible = 'off';
end
ha = axes();
% signalVar = zeros(size(img1Crop,3),2);
% noiseVar = zeros(size(img1Crop,3),2);
thrE_LTS = zeros(2,4, size(imgSCrop, 3));
fprintf('- Working on slice ')
strLoop = [];
for iS = 1:size(imgSCrop,3)
    %% feedback about analysed slice number on screen
    fprintf(repmat('\b',1,length(strLoop)));
    strLoop = sprintf('%d out of %d', iS, size(imgSCrop,3));
    fprintf('%s', strLoop);
    %% work with a single slice and all color channels
    C = squeeze(imgSCrop(:,:,iS,:));
    % C(:,:,3) = 0;

    %% calculate edges
    CE = C;
    method = 'Sobel';
    CE(:,:,1) = edge(C(:,:,1),method);
    CE(:,:,2) = edge(C(:,:,2),method);
        
    %% define thresholds using histogram for edge intensities
    %--- Extract intensities on edges froeach channel
    Red = C(:,:,1);
    Green = C(:,:,2);
    Red = Red(CE(:,:,1)==1);
    Green = Green(CE(:,:,2)==1);
    % Blue = C(:,:,3);
    %--- Get histValues for each channel
    [yRed, xR] = imhist(Red,128);
    [yGreen, xG] = imhist(Green,128);
    if ~all(xR==xG)
        error('Histogram binning for red and green is not matched!')
    end
    x = xR;
    % [yBlue, xB] = imhist(Blue);
    %--- Plot them together in one plot
    hold off
    plot(x, yRed, 'Red', x, yGreen, 'Green') %, xB, yBlue, 'Blue');
    hold on
    xlabel('intensity (normalized)')
    ylabel('count')
    ha.FontSize = 14;
    %--- fit Kernel distribution
    pdR = fitdist(Red,'Kernel','Support','positive');
    plot(x,pdR.pdf(x) * sum(yRed)/127 ,'k','LineWidth',1)
    pdG = fitdist(Green,'Kernel','Support','positive');
    plot(x,pdG.pdf(x) * sum(yGreen)/127 ,'k','LineWidth',1)
    
    thrE = zeros(2,4);
    
    %--- use percentiles as threshold levels
    x = 0:0.001:1;
    cdfT = pdR.cdf(x);
    %- 5th percentile
    [~, idx] = min((cdfT-5/100).^2);
    locL1 = x(idx);
    %- 20th percentile
    [~, idx] = min((cdfT-20/100).^2);
    locM = x(idx);
    %-----
    thrE(1,:) = [locL1 locM NaN NaN];
    
 
    %--- use percentiles as threshold levels
    cdfT = pdG.cdf(x);
    %- 25th percentile
    [~, idx] = min((cdfT-25/100).^2);
    locM = x(idx);
    %- 55th percentile
    [~, idx] = min((cdfT-55/100).^2);
    locU1 = x(idx);
    %-----
    thrE(2,:) = [NaN locM locU1 NaN];
    
    %--- plot threshold levels
    color='rg';
    for i=1:size(thrE,1)
        for j=1:size(thrE,2) -1
            plot([thrE(i,j) thrE(i,j)],ylim,color(i),'LineWidth',2)
            % text(thrE(i,j),diff(ylim)*(i/(3+j*0.1)),num2str(thrE(i,j)),'FontSize',14)
        end
    end
    %--- add text at the end to avoid obstruction by lines
    for i=1:size(thrE,1)
        for j=1:size(thrE,2) -1
            text(thrE(i,j),diff(ylim)*(i/(3+j*0.2)),num2str(thrE(i,j)),'FontSize',14)
        end
    end
    
    %---
    drawnow
    
    %% collect thresholds for all slices
    thrE_LTS(:,:,iS) = thrE;
    
end
fprintf('  ')
tee(fnameLog, display_time_delay)

close(hf)

    

%% smooth the thresholds
thrE_TLS_s = nan(size(thrE_LTS));
%x = (1:size(thrE_LTS,3))';
for iL = 1:size(thrE_LTS,1)
    for iT = 1:size(thrE_LTS,2)
        if ~all(isnan(thrE_LTS(iL,iT,:)))
            %--- smooth the curve
            y =squeeze(thrE_LTS(iL,iT,:));
            %--- smooth with moving average and find outliers
            thrE_TLS_s(iL,iT,:) = smooth_moving_average(y, 1);
        end
    end
end

%% show the smoothed thresholds
h = figure;
if figureHide
    h.Visible = 'off';
end
h.Position = [1032 537 752 808];
x = idxSliceOk(1):idxSliceOk(2);
subplot(2,1,1); hold on
plot(x, squeeze(thrE_LTS(1,:,:))','LineWidth',1)
plot(x, squeeze(thrE_TLS_s(1,:,:))','k')
set(gca,'FontSize',14)
title('Thresholds for red layer')
ylabel('intensity (normalized)')
subplot(2,1,2); hold on
plot(x, squeeze(thrE_LTS(2,:,:))','LineWidth',1)
plot(x, squeeze(thrE_TLS_s(2,:,:))','k')
set(gca,'FontSize',14)
title('Thresholds for green layer')
xlabel('slice')
ylabel('intensity (normalized)')
% legend({'lower 2/3 max','maximum','upper 2/3 max','upper 1/3 max'})
legend({'thr1','thr2','thr3','thr4'})

%% save picture
F = getframe(h);
imwrite(F.cdata, [trunk,'_prepro4_thresholds.png'])


%% Threshold the image
tee(fnameLog, 'Thresholding image slice by slice... \n')
tic
%--- Declare variable for thresholded 3D output volume
imgThr = zeros(size(imgSCrop),'uint8');
%---
fprintf('- Working on slice ')
strLoop = [];
for iS = 1:size(imgSCrop,3)
    %    clearvars -except iSlice img1Crop imgT img hdr imgSegmentation
    %     close all
    %% feedback about analysed slice number on screen
    fprintf(repmat('\b',1,length(strLoop)));
    strLoop = sprintf('%d out of %d', iS, size(imgSCrop,3));
    fprintf('%s', strLoop);
    
    %% work with a single slice and all color channels
    C = squeeze(imgSCrop(:,:,iS,:));
    thrE = thrE_TLS_s(:,:,iS);
    %% segment microglia cells
    %--- threshold the image
    BW1 = C(:,:,2) > thrE(2,3); %--- get green with a high threshold
    BW2 = C(:,:,2) > thrE(2,2) & C(:,:,1) < thrE(1,1); %--- green with low threshold (maximum of int at edges) but excluding red
    BW = BW1 | BW2;
       
    %% segment nuclei
    BWnuclei = C(:,:,1) > thrE(1,2) & C(:,:,2) > thrE(2,3); %--- get the nuclei
    
    %% insert thresholded slice into 3D volume
    imgThr(:,:,iS, 1) = BWnuclei;
    imgThr(:,:,iS, 2) = BW;
    
    
end
fprintf('  ')
tee(fnameLog, display_time_delay)

%% remove small cluster from green layer
tee(fnameLog, 'Removing small cluster from green layer ...  ')
tic
%--- clusterize
L = bwlabeln(imgThr(:,:,:,2));
%--- get area of cluster
T = regionprops('table', L, 'centroid', 'area');
%--- x/y dimension are exchanged by "regionprops"
S = [];
S.cogI = T.Centroid(:,[2 1 3]);
S.cog = S.cogI .* repmat(hdr.pixdim(2:4), height(T), 1) - repmat(hdr.pixdim(2:4)./2, height(T), 1);
S.voxN = T.Area;
S.voxV = T.Area * prod(hdr.pixdim(2:4));
%--- find small cluster with threshold
S.smallCluster = S.voxV < thrVolume;
%--- remove small cluster
idxBig = find(S.smallCluster==0);
IG = ismember(L,idxBig);
tee(fnameLog, display_time_delay)


%% %% Fill holes
tee(fnameLog, 'Fill holes ingreen layer ...  ')
tic
IG = imfill(IG,'holes');
% IGF = imfill(IG,'holes');
tee(fnameLog, display_time_delay)


%% Define the soma of the cells, by removing the branches from it
%--- create structuring element
%- radius of sphere
rSomaPreliminary = 1.5;
dim = hdr.pixdim(2:4);
sesz = rSomaPreliminary ./ dim;
%--- center of the sphere
secr = round(sesz) + 1;
%--- size of bounding box
sesz = 2 * round(sesz) + 1;
%--- distance from center
ne = zeros(sesz);
for i=1:size(ne,1)
    for j=1:size(ne,2)
        for k=1:size(ne,3)
            ne(i,j,k) = sum(([i j k] .* dim - secr .* dim).^2).^.5;
        end
    end
end
%--- voxel within requested radius from center
ne(ne>rSomaPreliminary) = 0;
ne(secr(1),secr(2),secr(3)) = 1;
ne = logical(ne);
%--- visualize the spherical structuring element
% show_cells_3D(ne,hdr,[],[],10);
%---
seSoma = strel(ne);
%--- open the green layer to remove arms
tee(fnameLog, 'Define soma  ...  ')
tic
IS = imopen(IG,seSoma);
tee(fnameLog, display_time_delay)

%% Fill small gaps in MaskCell (green layer) in the immediate vicinity of the soma
%- Two-step procedure:
%- 1. Dilate the soma (smoothed by imopen above) to restrict the gap filling operation to the immediate surrounding of the soma
%- 2. Fill gaps in MaskCell, using the Matlab function imclose with a spherical structuring element of radiaus 1.2 microns
%--- Create structuring element
rDilate = 1.2;
dim = hdr.pixdim(2:4);
sesz = rDilate ./ dim;
%--- center of the sphere
secr = round(sesz) + 1;
%--- size of bounding box
sesz = 2 * round(sesz) + 1;
%
ne = zeros(sesz);
for i=1:size(ne,1)
    for j=1:size(ne,2)
        for k=1:size(ne,3)
            ne(i,j,k) = sum(([i j k] .* dim - secr .* dim).^2).^.5;
        end
    end
end
ne(ne>rDilate) = 0;
ne(secr(1),secr(2),secr(3)) = 1;
ne = logical(ne);
%---
% show_cells_3D(ne,hdr,[],[],10);
%---
seDilate = strel(ne);
%--- dilate the soma
tee(fnameLog, 'Dilate soma mask ...  ')
tic
ID = imdilate(IS,seDilate);
tee(fnameLog, display_time_delay)
%--- imclose the green layer
tee(fnameLog, 'Close small gaps in green layer ...  ')
tic
IGc = imclose(IG,seDilate);
tee(fnameLog, display_time_delay)
%--- mask closed green layer with dilated soma to restrict the closing operation to the surrounding of the soma
ID = IGc & ID;
%--- add the closed soma to the green layer
IG = IG | ID;


%% Mask the nuclei exclusively with the MaskSoma to exclude voxel cluster which might be artefacts. 
%--- For example, if a branch lies in close proximity to a non-microglia nucleus, 
%--- an overlap of red and green might result, erroneously indicating the presence
%--- of a microglia nucleus in MaskNuclei.
IN = imgThr(:,:,:,1) & IS;

%% remove small cluster of potential nuclei to exclude those, which are actually not a soma but have arteficial red pixel in it
%--- clusterize
L = bwlabeln(IN);
%--- get area of cluster
T = regionprops('table', L, 'centroid', 'area');
%--- x/y dimension are exchanged by "regionprops"
S = [];
S.cogI = T.Centroid(:,[2 1 3]);
S.cog = S.cogI .* repmat(hdr.pixdim(2:4), height(T), 1) - repmat(hdr.pixdim(2:4)./2, height(T), 1);
S.voxN = T.Area;
S.voxV = T.Area * prod(hdr.pixdim(2:4));
%--- find small cluster with threshold
S.smallCluster = S.voxV < thrVolume;
%--- remove small cluster
idxBig = find(S.smallCluster==0);
IN = ismember(L,idxBig);

%% Now redefine the soma in the repaired mask, by removing the branches from it, but using a larger diameter
%--- create structuring element
%- radius of sphere
dim = hdr.pixdim(2:4);
sesz = rSoma ./ dim;
%--- center of the sphere
secr = round(sesz) + 1;
%--- size of bounding box
sesz = 2 * round(sesz) + 1;
%--- distance from center
ne = zeros(sesz);
for i=1:size(ne,1)
    for j=1:size(ne,2)
        for k=1:size(ne,3)
            ne(i,j,k) = sum(([i j k] .* dim - secr .* dim).^2).^.5;
        end
    end
end
%--- voxel within requested radius from center
ne(ne>rSoma) = 0;
ne(secr(1),secr(2),secr(3)) = 1;
ne = logical(ne);
%--- visualize the spherical structuring element
% show_cells_3D(ne,hdr,[],[],10);
%---
seSoma = strel(ne);
%--- open the green layer to remove arms
tee(fnameLog, 'Redefine soma after gap closure ...  ')
tic
IS = imopen(IG,seSoma);
tee(fnameLog, display_time_delay)


%% Add nuclei to the soma mask
%- to account for the possibility that a nucleus is not anymore inside a soma, after redifining them with a different diameter
IS = IS | IN;

%% accept soma only, if nucleus can be found within
L = bwlabeln(IS);
idxSoma = unique(L(IS & IN));
IS = ismember(L,idxSoma);


%% Dilate the redefined soma, for following purposes:
%- 2. separate real branches from minor bumps on its surface
%- 3. separate branches which share the same basis or are connected by ridges on the surface of the soma
%---
tee(fnameLog, 'Dilate soma ...  ')
tic
ID = imdilate(IS,seDilate);
ID = IG & ID;
tee(fnameLog, display_time_delay)

%% accept dilated-soma only, if nucleus (soma) can be found within (this step is needed, because during the dilation step, the dilation can swap across background and create disconnected areas)
L = bwlabeln(ID);
idxSomaDil = unique(L(ID & IN));
ID = ismember(L,idxSomaDil);

%% final result of segmentation: arms + border + soma + nuclei
%--- summing up results in following labels:
%- arms   = 1
%- border = 2
%- soma   = 3
%- nuclei = 4
IABSN = uint8(IG + ID + IS + IN);
cmap = [0 0 0; 0 0.8 0; 0.8 0.3 1; 0 0 1; 1 1 0];
%--- show final result
h = show_cells_3D(IABSN,hdr,[],[],figureScaling, cmap, figureHide);

%--- save picture
F = getframe(h);
imwrite(F.cdata, [trunk,'_ortho4_segmented.png']);

%% histogram of intensities within the nuclei and microglia segments
%--- create mask
if idxSliceOk(1) > 1
    temp1 = zeros(size(IN,1), size(IN,2), idxSliceOk(1)-1, 2);
else
    temp1 = [];
end
if idxSliceOk(2) < size(img1,3)
    temp2 = zeros(size(IN,1), size(IN,2), size(img1,3)-idxSliceOk(2), 2);
else
    temp2 = [];
end
mask = cat(3, temp1, cat(4, IN, IG), temp2);
%---
hf = volume_histogram_by_slice_and_color(img1, 20, mask, [], [], figureHide);
%% save picture
%--- histogram
F = getframe(hf(1));
imwrite(F.cdata, [trunk,'_prepro5_histograms_segmented.png'])

%%
tee(fnameLog, 'Save images ...  ')
tic
%% save data: equalized image
fnameOut = [trunk, '_stack2_equalized.mat'];
%--- convert image to 8 bit color depth
write_3d_image_matrix(hdr, uint8(imgEq*255), fnameOut)

%% save data: final segmentation
fnameOut = [trunk, '_stack3_segmented.mat'];
hdr2 = change_hdr_datatype_info(hdr, 2);
hdr2 = change_hdr_dim(hdr2, size(IABSN));
hdr2.cmap = cmap;
write_3d_image_matrix(hdr2, uint8(IABSN), fnameOut)

%%
tee(fnameLog, display_time_delay)

%% Make a video
hf = figure;
if figureHide
    hf.Visible = 'off';
end
r = groot;
tee(fnameLog, 'Make movie ...')
tic
movName = [trunk,'_video_soma_cell'];
v = VideoWriter(movName, 'MPEG-4');
v.FrameRate = 9;
open(v)
for iS=1:size(IABSN,3)
    temp = rot90(squeeze(IABSN(:,:,iS)));
    temp = ind2rgb(temp,cmap);
    C = rot90(squeeze(imgEqCrop(:,:,iS,:)));
    %--- normalize to a certain percentile range
    prct = prctile(reshape(C(:,:,1),[],1),[1 99]);
    C(:,:,1) = C(:,:,1) - prct(1);
    C(:,:,1) = C(:,:,1)./diff(prct);
    prct = prctile(reshape(C(:,:,2),[],1),[1 99]);
    C(:,:,2) = C(:,:,2) - prct(1);
    C(:,:,2) = C(:,:,2)./diff(prct);
    C(C<0) = 0;
    C(C>1) = 1;
    %---
    C(:,:,3) = 0;
    %--- flip colors
    C = C(:,:,[3 2 1]);
    %---
    temp3 = ones(size(IABSN,1), 2, 3);
    r.CurrentFigure = hf;
    imshow(cat(2, C, temp3, temp), [])
    writeVideo(v,getframe(hf))
end
close(v)
close(hf)
tee(fnameLog, display_time_delay)

%%
tee(fnameLog, 'Segmentation completed!\n')
tee(fnameLog, display_time_delay(hTimer))
fprintf('\n\n')
