function [h, D] = show_cells_3D(img, hdr, sliceN, clPrct, scaleFactor, cmap, hide)
% function [h, D] = show_cells_3D(img, hdr, sliceN, clPrct, scaleFactor, cmap, hide)
% function [h, D] = show_cells_3D(fnameImgage, [], sliceN, clPrct, scaleFactor, cmap, hide)
% Input arguments:
%   img         - 3D or 4D image matrix (4. dimension corresponds to color channel; max. 3 channels); "img" is sufficient as input.
%   fnameImage  - Path to an image file; "fnameImage" is sufficient as input
%   hdr         - Header conforming to the format of a nifti-1 file header (optional; but important for correct scaling of image dimensions)
%   sliceN      - 1x3 vector with the cursor coordinates to show after opening the figure (optional)
%   clPrct      - 1x2 vector with the lower and upper percentiles for intensity clipping (optional)
%   scaleFactor - scalar for size scaling of the image, (optional; by default it
%                 has the value 1 and means that each pixel on the screen corresponds
%                 to one unit (e.g. micrometer) as calculated using the scaling info in "hdr"
% Output arguments:
%   h           - handle to the created figure
%   D           - structure with fields for all data used during visualization
%
% Please note: 
% For speed during visualization, the bit depth of the "img" data-type is crucial.
% Therefore RGB images will be converted always into indexed images with
% colormap. If you provide grey scale images, it is best practice to provide
% an integer data-type with the lowest possible bit depth (e.g. uint8).
%
%% Announce me in the command window
fprintf('Creating view of orthogonal slices:\n');
%--- check whether the figure should be hidden (this can be useful, if only a screenshot of the figure is desired)
if exist('hide', 'var') && (islogical(hide) || isnumeric(hide)) && hide
    fprintf('-Figure will be hidden\n')
    fprintf('-Note: Use the returned handle to make this figure visible or to close it\n')
end
%% default  parameters
%--- color order
colorOrder = [3 2 1]; %--- make first layer blue, second green, third red
% colorOrder = [2 3 1]; %--- make first layer blue, second red, third green
% colorOrder = [2 1 3]; %--- make first layer green, second red, third blue
%--- convert RGB to color map 
%- This option significantly increases speed of visualization, when moving the cursor over the axes. 
%- However, it can take quite long for the conversion step at the beginning.
%- Also, precision will be lost. If reading out exact intensity values is important to you, don't convert to map
% convert2cmap = false;
convert2cmap = true;
%% delete empty input variables
inputVars = {'img','hdr','sliceN','clPrct', 'scaleFactor', 'cmap', 'hide'};
for i = 1:nargin
    if isempty(eval(inputVars{i}))
        eval(['clearvars ', inputVars{i}])
    end
end

%% check whether to yoke to other figures is true
D.yoke = true; %--- this paramter is fixed at the moment; consider to insert it as input argument

%% check whether "img" is a file name, in case load image from file
if ischar(img)
    [fname, ext] = regexp(img,'\.nii(\.gz)*$|\.mat$','split','match');
    if isempty(ext) || ~ismember(ext{1}, {'.nii','.nii.gz','.mat'})
        error('Image-Filename has no valid extension. Allowed are: ''.mat|.nii|.nii.gz''')
    else
        fprintf('-Reading image-file:\n   %s\n', img)
        tic
        switch ext{1}
            case {'.nii','.nii.gz'}
                if ~exist([fname{1} '.nii'], 'file') && ~exist([fname{1} '.nii.gz'], 'file')
                    error('Neither compressed nor uncompressed nifti file found: %s.nii(.gz)', fname{1})
                else
                    [hdr, img] = read_nifti_gz(img);
                end
            case '.mat'
                if ~exist(img, 'file')
                    error('File not found: %s', img)
                else
                    temp = load(img);
                    %--- check whether the file contained a structure
                    tempFN = fieldnames(temp);
                    if length(tempFN)==1 && isstruct(temp.(tempFN{1}))
                        temp = temp.(tempFN{1});
                    end
                    if isfield(temp,'img')
                        img = temp.img;
                    else
                        error('Provided ''.mat'' file has to contain a variable ''img''')
                    end
                    if isfield(temp,'hdr')
                        hdr = temp.hdr;
                    end
                    if isfield(temp,'cmap')
                        cmap = temp.cmap;
                    end
                end
        end
        toc
    end
elseif islogical(img) &&  length(size(img))>=3
    img = uint8(img);
elseif ~isnumeric(img) || length(size(img))<3
    error('First argument has to be an 3D image matrix or valid image file-path (e.g.: .nii/.nii.gz/.mat)')

end

%% check if cmap exists either as variable or field of the header
if ~exist('cmap', 'var')
    if exist('hdr','var') && isstruct(hdr) && isfield(hdr,'cmap')
        cmap = hdr.cmap;
    end
end


%% define indices of slices to be shown
sz = size(img);
if ~exist('sliceN','var')
    D.sliceN = round(sz(1:3)/2);
elseif length(sliceN)==1
    D.sliceN = repmat(sliceN,3,1);
else
    D.sliceN = sliceN;    
end
idxOutOfRange = D.sliceN > sz(1:3);
if any(idxOutOfRange)
    warning('%s       \n%s', 'The selected cross-hair position is out of range!', 'Using the largest possible index instead!')
    sliceNt = sz(1:3);
    D.sliceN(idxOutOfRange) = sliceNt(idxOutOfRange);
end
%% Check whether image has any content
if range(img(:))==0
    error(sprintf('\nRange of values in image equals 0!\n => Nothing to visualize!!!\n    Aborting program!!!'))
end

%% Apply intensity scaling, depending on image dimensionality (4D image will be interpreted as RGB)
%--- 3D means grey scale image or indexed image
%--- 4D means color information given, 4th dimension corresponds to color channel
if length(sz)==3 && exist('cmap','var') %--- 3D indexed image with color map
    fprintf('-Image was found to be indexed with a color map\n')
    cl = [0 size(cmap,1)-1];
    D.marker = '+r';    
elseif length(sz)==3 && exist('clPrct', 'var') %--- 3D gray scale
    fprintf('-Clipping data as requested [%g %g]...',clPrct)
    tic
    %--- clipping with function "imagesc" works only for 3D but not for 4D images
    if sum(img(:)==0) > numel(img)*(5/100) %--- if proportion of zeros > 5%
        cl = prctile(double(img(img~=0)),clPrct);
    else
        cl = prctile(double(img(:)),clPrct);
    end
    if cl(1)==cl(2)
        cl(1) = min(img(:));
        cl(2) = max(img(:));
    end
    toc
    D.marker = '+r';
elseif length(sz)==3
    clPrct = [0 100];
    cl(1) = min(img(:));
    cl(2) = max(img(:));
    D.marker = '+r';
elseif length(sz)==4 %--- 4D
    img = double(img); %--- convert to double
    if exist('clPrct', 'var') %--- clip only if requested
        fprintf('-Clipping data as requested [%g %g]...',clPrct)
        tic
        for i=1:sz(4)
            
            %--- clip each color layer at desired percentiles and scale separately to yield max==1
            temp = img(:,:,:,i);
            pcl = prctile(temp(:),clPrct); %--- prctile works only with double precision
            temp(temp>pcl(2)) = pcl(2);
            temp(temp<pcl(1)) = pcl(1);
            temp = temp - pcl(1);
            img(:,:,:,i) = temp ./ max(temp(:));
            %--- scale each color layer separately to yield max==1 (without clipping)
            % img(:,:,:,i) = img(:,:,:,i) / max(max(max(img(:,:,:,i))));
        end
        toc
    else %--- allowed intensity range for rgb images is from 0 to 1
        imgMin = min(img(:));
        if imgMin < 0
            img = img - imgMin;
            D.colorScaling = imgMin;
        else
            D.colorScaling = 0;
        end
        imgMax = max(img(:));
        if imgMax >1
            img = img ./ imgMax;
            D.colorScaling(2) = imgMax;
        else
            D.colorScaling(2) = 1;
        end
        if all(D.colorScaling == [0 1])
            D = rmfield(D,'colorScaling');
        end
    end
    if sz(4)==2
        img(:,:,:,3) = 0;
    elseif sz(4)>3
        fprintf('-WARNING: Your image has %d color channels!\n', sz(4));
        fprintf('          Only the first three channels can be shown as RGB!\n')
        img = img(:,:,:,1:3);
    end
    cl = [0.1 0.2];
    D.marker = '+w';
else
    errorStr = repmat('%d ',1,length(sz));
    errorStr = sprintf('Unexpected dimensionality of image: %s\n',errorStr);
    error(errorStr, sz);
end

%% Adjust color order
if exist('colorOrder','var') && ~isempty(colorOrder) && size(img,4)==3
    img = img(:,:,:,colorOrder);
    
end
    
%% Converting RGB to colormap
if exist('convert2cmap','var') && convert2cmap && size(img,4)==3
    fprintf('-Converting RGB to colormap...')
    tic
    imgT = reshape(img, size(img,1)*size(img,2),size(img,3),size(img,4));
    if ~any(strcmpi(class(imgT),{'single','double'})) %--- rgb2ind doesn't work with integers
        imgT = double(imgT);
    end
    [imgT,cmap] = rgb2ind(imgT,2^16);
    % cmap = rgb2gray(cmap); %--- to gray scale
    img = reshape(imgT, size(img,1),size(img,2),size(img,3));
    cl = [0 size(cmap,1)-1];
    toc
end

%% check/optimize data type of inidces to color map
if exist('cmap','var')  
    Lmax = max(img(:));
    typeStr = {'uint8','uint16','uint32','uint64'};
    for i=1:length(typeStr)
        if strcmp(class(img),typeStr{i})
            break
        elseif Lmax <= intmax(typeStr{i})
            fprintf('-Optimizing data type for faster visualization: ''%s'' to ''%s''\n', class(img), typeStr{i})
            img = eval([typeStr{i},'(img)']);
            break
        end
    end
end

%% get screen dimension
innerPosition = get( groot, 'Screensize' );
%--- note: note the whole screen can be used, therefore limit the useable part of the screen
innerPosition([3 4]) = innerPosition([3 4])-[60, 110];



%% get image and pixel dimension (scaling)
if ~exist('hdr','var')
    D.pix = ones(1,3);
elseif ~isstruct(hdr) && isnumeric(hdr) && isvector(hdr) && length(hdr)<=3
    D.pix = ones(1,3);
    D.pix(1:length(hdr)) = hdr;
elseif isstruct(hdr) && isfield(hdr, 'pixdim')
    D.pix = hdr.pixdim(2:4);
end
szImg = size(img);
if length(szImg) < 3 %--- stacks with only 1 slice (2D image)
    xyz = [szImg 1] .* D.pix;
else
    xyz = szImg(1:3) .* D.pix;
end
D.pixH = D.pix/2;

%% define figure size 
border = 15; %--- border between figure next and between axes, in pixel
bw = xyz(1)+xyz(3)+ 3*border; %--- breite window
hw = xyz(2)+xyz(3)+ 3*border; %--- height window
bwr = bw/(innerPosition(3)); %--- ratio breite window to breite screen
hwr = hw/(innerPosition(4)); %--- ratio height window to hight screen

%% scale figure if requested
if exist('scaleFactor', 'var') 
    if isnumeric(scaleFactor)
        fprintf('-Scaling image size by user specified factor %g\n', scaleFactor)
        bw = bw * scaleFactor;
        hw = hw * scaleFactor;
        bwr = bw/(innerPosition(3)); %--- ratio breite window to breite screen
        hwr = hw/(innerPosition(4)); %--- ratio height window to hight screen
    elseif ischar(scaleFactor) && strcmpi(scaleFactor, 'optimize')
        %--- if the ratio is bigger than 1, then scale the window
        if bwr<0.90 && hwr<0.90
            if bwr>hwr
                bw = bw/bwr;
                hw = hw/bwr;
                fprintf('-NOTE: Image was scaled by %.2f to optimize screen usage\n', 1/bwr)
            else
                bw = bw/hwr;
                hw = hw/hwr;
                fprintf('-NOTE: Image was scaled by %.2f to optimize screen usage\n', 1/hwr)
            end
            bwr = bw/(innerPosition(3)); %--- ratio breite window to breite screen
            hwr = hw/(innerPosition(4)); %--- ratio height window to hight screen
        end
    else
        fprintf('-Scale factor: ')
        disp(scaleFactor)
        error('Forbidden scaling parameter chosen (see above)')
    end
end

%% check whether figure size fits on the screen, otherwise scale it down
if bwr>1.001 || hwr>1.001
    if bwr>hwr
        bw = bw/bwr;
        hw = hw/bwr;
        fprintf('-WARNING: Image was scaled by %.2f to fit on screen\n', 1/bwr)
    else
        bw = bw/hwr;
        hw = hw/hwr;
        fprintf('-WARNING: Image was scaled by %.2f to fit on screen\n', 1/hwr)
    end
    
end
    
%% define 
%--- subtract borders from window size, to get the drawable part of the window
bwc = bw - 3*border;
hwc = hw - 3*border;

b1 = xyz(1)/(xyz(1)+xyz(3)) * bwc / bw; %--- Breite of left panel
b2 = xyz(3)/(xyz(1)+xyz(3)) * bwc / bw; %--- Breite of right panel
h1 = xyz(2)/(xyz(2)+xyz(3)) * hwc / hw; %--- Height of lower panel
h2 = xyz(3)/(xyz(2)+xyz(3)) * hwc / hw; %--- Height of upper panel

bb = (1-(b1+b2))/3; %--- border in x direction
hb = (1-(h1+h2))/3; %--- border in y direction

%--- XZ (oben links)
positionVector(1,:) = [bb hb*2+h1 b1 h2]; %--- links, unten, breite, hoehe
%--- ZY (unten rechts)
% positionVector(2,:) = [0.1+c1b 0.1+r2h c2b r1h]; %--- links, unten, breite, hoehe
positionVector(2,:) = [2*bb+b1 hb b2 h1]; %--- links, unten, breite, hoehe
%--- XY (unten links)
positionVector(3,:) = [bb hb b1 h1]; %--- links, unten, breite, hoehe

%%
h = figure('Units','Pixel','Position',[60 10 bw hw]);
if exist('hide', 'var') && (islogical(hide) || isnumeric(hide)) && hide
    h.Visible = 'off';
end
h.Pointer = 'crosshair';

%% define user interface and callbacks
set (h, 'WindowButtonMotionFcn', @mouseMove);
set (h, 'WindowButtonDownFcn', @mouseDown);
set (h, 'WindowButtonUpFcn', @mouseUp);
%--- text box to edit cursor position
hpos = uicontrol('Style', 'edit',...
    'Position', [20 0 130 15],...
    'Callback', @text_edit,...
    'TooltipString', 'coordinates (voxels)');
hpos.String = sprintf('%d %d %d',D.sliceN);
%--- text box to edit percentile range for clipping
if exist('clPrct', 'var')
hClPrct = uicontrol('Style', 'edit',... 
    'Position', [160 0 60 15],...
    'Callback', @text_edit_clPrct,...
    'TooltipString', 'percentile range');
    hClPrct.String = sprintf('%d %d',clPrct);
end
%%
%--- define a title
% if isfield(hdr,'datatypestr')
%     titleStr = hdr.datatypestr;
%     set(h, 'Name', titleStr, 'NumberTitle','off')
% end

if isfield(hdr,'name')
    D.titleStr = [hdr.name ';'];
else
    D.titleStr = sprintf('Dimension: %d %d %d;', sz(1:3));
end
h.Name = sprintf('%s Point: %d %d %d',D.titleStr, D.sliceN);
h.NumberTitle = 'off';


%% positioning of images within axes coorinates
Xpos = [1 sz(1)]*D.pix(1)-D.pixH(1);
Ypos = [1 sz(2)]*D.pix(2)-D.pixH(2);
Zpos = [1 sz(3)]*D.pix(3)-D.pixH(3);


%% plot the orthogonal image planes
%p1 = subplot(2,2,1)
axs(1) = subplot('Position',positionVector(1,:));
hold on
D.hI(1) = imagesc('XData', Xpos,'YData', Zpos, 'CData',permute(squeeze(img(:,D.sliceN(2),:,:)),[2 1 3]),cl); %--- X/Z
axis tight
axis equal
D.hL(1) = plot(      repmat(D.sliceN(1)*D.pix(1)-D.pixH(1),2,1), ylim, D.marker(end));
D.hL(2) = plot(xlim, repmat(D.sliceN(3)*D.pix(3)-D.pixH(3),2,1),       D.marker(end));
%--- add scaling markers corresponding to 10 units
% plot([1, 10],[1 1],'r','LineWidth',3)
% plot([1 1],[1, 10],'r','LineWidth',3)
% axis off
% p2 = subplot(2,2,2)
axs(2,1) = subplot('Position',positionVector(2,:));
hold on
D.hI(2) = imagesc('XData', Zpos,'YData', Ypos, 'CData',squeeze(img(D.sliceN(1),:,:,:)),cl); %--- Z/Y
axis tight
axis equal
D.hL(3) = plot(      repmat(D.sliceN(3)*D.pix(3)-D.pixH(3),2,1), ylim, D.marker(end));
D.hL(4) = plot(xlim, repmat(D.sliceN(2)*D.pix(2)-D.pixH(2),2,1),       D.marker(end));
%--- add scaling markers corresponding to 10 units
% plot([1, 10],[1 1],'r','LineWidth',3)
% plot([1 1],[1, 10],'r','LineWidth',3)
% axis off
% p3 = subplot(2,2,3)
axs(3,1) = subplot('Position',positionVector(3,:));
hold on
D.hI(3) = imagesc('XData', Xpos,'YData', Ypos, 'CData', permute(squeeze(img(:,:,D.sliceN(3),:)),[2 1 3]),cl); %--- X/Y
axis tight
axis equal
D.hL(5) = plot(      repmat(D.sliceN(1)*D.pix(1)-D.pixH(1),2,1), ylim, D.marker(end));
D.hL(6) = plot(xlim, repmat(D.sliceN(2)*D.pix(2)-D.pixH(2),2,1),       D.marker(end));
%--- add scaling markers corresponding to 10 units
% plot([1, 10],[1 1],'r','LineWidth',3)
% plot([1 1],[1, 10],'r','LineWidth',3)
% axis off
% axis equal
h.Pointer = 'crosshair';
%---
if exist('cmap','var')
    colormap(cmap);
else
    colormap(gray(256))
end
drawnow
%%
D.axs = axs;
D.img = img;
D.hdr = hdr;
if exist('cmap','var')
D.cmap = cmap;
end
D.positionVector = positionVector;
D.mouseState = false;
D.cl = cl;
if ~exist('clPrct','var')
    clPrct = [1 100];
end
D.clPrct = clPrct;
D.hpos = hpos;
h.UserData = D;
%%
% figure;
% [X,Y]=meshgrid(1:20,1:20);
% image(Y)

%%
if nargout==2
    fields2keep = {'img','hdr','cmap'};
    fields2remove   = fieldnames(D);
    fields2remove(ismember(fields2remove, fields2keep)) = [];
    for i=1:length(fields2remove)
        D = rmfield(D,fields2remove{i});
    end
end


%% end of function

end
%%
function mouseMove(hObject, eventdata)
D = hObject.UserData;
if D.mouseState==true
idx = 0;
for i=1:3
    C = get(D.axs(i), 'CurrentPoint');
    % C = round(C(1,1:2));
    C = C(1,1:2);
    
    xl = xlim(D.axs(i));
    yl = ylim(D.axs(i));
    if C(1) >= xl(1) && C(1) <= xl(2) && C(2) >= yl(1) && C(2) <= yl(2)
        idx = i;
        break
    end
end
if idx>0
    switch idx
        case 1
            %D.sliceN([1 3]) = C;
            D.sliceN([1 3]) = round((C + D.pixH([1 3])) ./D.pix([1 3]));
        case 2
            % D.sliceN([3 2]) = C;
            D.sliceN([3 2]) = round((C + D.pixH([3 2])) ./D.pix([3 2]));
        case 3
            % D.sliceN([1 2]) = C;
            D.sliceN([1 2]) = round((C + D.pixH([1 2])) ./D.pix([1 2]));
    end
    hObject.UserData = D;
    D.hpos.String = sprintf('%d %d %d', D.sliceN);
    update_views(hObject)
end
end
end

function mouseDown(hObject, eventdata)
hObject.UserData.mouseState = true;
mouseMove(hObject, [])
end

function mouseUp(hObject, eventdata)
D = hObject.UserData;
D.hpos.String = sprintf('%d %d %d',D.sliceN);
hObject.UserData.mouseState = false;
end

function text_edit(hObject, eventdata, foreignAction)
xyzStr = get(hObject,'String');
xyz = cellfun(@str2double,regexp(xyzStr, ' ', 'split'));
if length(xyz)<3 || any(isnan(xyz))
    hObject.String = 'Three integers needed!';
else
    h = hObject.Parent;
    sz = size(h.UserData.img);
    xyz(xyz<1) = 1;
    xyz(xyz>sz(1:3)) = sz(xyz>sz(1:3));
    set(hObject,'String',sprintf('%d %d %d', xyz));
    h.UserData.sliceN = xyz;
    if exist('foreignAction', 'var')
        update_views(h, foreignAction)
    else
        update_views(h)
    end
end
end


function text_edit_clPrct(hObject, eventdata)
h = hObject.Parent;
clPrctStr = get(hObject,'String');
clPrct = cellfun(@str2double,regexp(clPrctStr, ' ', 'split'));
sz = size(h.UserData.img);
if isfield(h.UserData,'cmap') && sz~=3
    errordlg({'Sorry, changing percentile range interactively is implemented only for single-layer grey scale images!'})
    hObject.String = sprintf('%d %d', h.UserData.clPrct);
end
if length(clPrct)<2 || any(isnan(clPrct))
    hObject.String = 'Two integers needed!';
else
    if h.UserData.clPrct == clPrct
        errordlg({'The values did not change!'})
    else
        %---
        fprintf('-Clipping data as requested [%g %g]...',clPrct)
        tic
        %--- calculate the percentiles across the whole 3D volume
        img = h.UserData.img;
        if sum(img(:)==0) > numel(img)*(5/100) %--- if proportion of zeros > 5%
            cl = prctile(double(img(img~=0)),clPrct);
        else
            cl = prctile(double(img(:)),clPrct);
        end
        %---
        % if cl(1)==cl(2)
        %     cl(1) = min(img(:));
        %     cl(2) = max(img(:));
        % end
        %---
        h.UserData.cl = cl;
        %--- assign new clipping values
        for i=1:length(h.UserData.axs)
            h.UserData.axs(i).CLim = cl;
        end
        %---
        % update_views(h)
        %---
        toc
    end
end
end

% function update_views(h)
% D = h.UserData;
% % fprintf('Point: %d %d %d\n', D.sliceN)
% %---
% axes(D.axs(1));
% imagesc(rot90(squeeze(D.img(:,D.sliceN(2),:,:))),D.cl) %--- X/Z
% % imagesc(rot90(squeeze(D.img(:,D.sliceN(2),:,:)))) %--- X/Z
% axis off
% hold on
% % plot(D.sliceN(1),size(D.img,3)+1-D.sliceN(3),D.marker, 'MarkerSize', 20)
% plot(repmat(D.sliceN(1),2,1),ylim,D.marker(end))
% plot(xlim,repmat(size(D.img,3)+1-D.sliceN(3),2,1),D.marker(end))
% hold off
% axes(D.axs(2));
% imagesc(flipud(squeeze(D.img(D.sliceN(1),:,:,:))),D.cl) %--- Y/Z
% % imagesc(flipud(squeeze(D.img(D.sliceN(1),:,:,:)))) %--- Y/Z
% axis off
% hold on
% % plot(D.sliceN(3),size(D.img,2)+1-D.sliceN(2),D.marker, 'MarkerSize', 20)
% plot(repmat(D.sliceN(3),2,1),ylim,D.marker(end))
% plot(xlim,repmat(size(D.img,2)+1-D.sliceN(2),2,1),D.marker(end))
% hold off
% % p3 = subplot(2,2,3)
% % D.axs(3,1) = subplot('Position',D.positionVector(3,:));
% axes(D.axs(3));
% imagesc(rot90(squeeze(D.img(:,:,D.sliceN(3),:))),D.cl) %--- X/Y
% % imagesc(rot90(squeeze(D.img(:,:,D.sliceN(3),:)))) %--- X/Y
% axis off
% hold on
% % plot(D.sliceN(1),size(D.img,2)+1-D.sliceN(2),D.marker, 'MarkerSize', 20)
% plot(repmat(D.sliceN(1),2,1),ylim,D.marker(end))
% plot(xlim,repmat(size(D.img,2)+1-D.sliceN(2),2,1),D.marker(end))
% hold off
% %---
% if isfield(D,'cmap')
% colormap(D.cmap);
% end

function update_views(h, foreignAction)
D = h.UserData;
% fprintf('Point: %d %d %d\n', D.sliceN)
%---
for i=1:length(D.hL)
    delete(D.hL(i))
end
%---
% axes(D.axs(1));
D.hI(1).CData = permute(squeeze(D.img(:,D.sliceN(2),:,:)),[2 1 3]); %--- X/Z
D.hL(1) = plot(D.axs(1),      repmat(D.sliceN(1)*D.pix(1)-D.pixH(1),2,1), ylim(D.axs(1)), D.marker(end));
D.hL(2) = plot(D.axs(1),xlim(D.axs(1)), repmat(D.sliceN(3)*D.pix(3)-D.pixH(3),2,1),       D.marker(end));

% axes(D.axs(2));
D.hI(2).CData = squeeze(D.img(D.sliceN(1),:,:,:)); %--- Z/Y
D.hL(3) = plot(D.axs(2),      repmat(D.sliceN(3)*D.pix(3)-D.pixH(3),2,1), ylim(D.axs(2)), D.marker(end));
D.hL(4) = plot(D.axs(2),xlim(D.axs(2)), repmat(D.sliceN(2)*D.pix(2)-D.pixH(2),2,1),       D.marker(end));

% axes(D.axs(3));
D.hI(3).CData = permute(squeeze(D.img(:,:,D.sliceN(3),:)),[2 1 3]); %--- X/Y
D.hL(5) = plot(D.axs(3),      repmat(D.sliceN(1)*D.pix(1)-D.pixH(1),2,1), ylim(D.axs(3)), D.marker(end));
D.hL(6) = plot(D.axs(3),xlim(D.axs(3)), repmat(D.sliceN(2)*D.pix(2)-D.pixH(2),2,1),       D.marker(end));

%---
%axes(D.axs(idx));
index = [];
intensity = D.img(D.sliceN(1),D.sliceN(2),D.sliceN(3),:);
if isfield(D,'cmap')
    index = intensity;
    intensity = D.cmap(intensity+1,:);
end
if isfield(D,'colorScaling')
    intensity = (intensity .* D.colorScaling(2)) + D.colorScaling(1);
end
if isempty(index)
    formatStr = ['%s Point: %d %d %d; Intensity:', repmat(' %.3g',1,length(intensity))];
    h.Name = sprintf(formatStr, D.titleStr, D.sliceN, intensity);
else
    formatStr = ['%s Point: %d %d %d; Index: %d; Intensity:', repmat(' %.3g',1,length(intensity))];
    h.Name = sprintf(formatStr, D.titleStr, D.sliceN, index, intensity);
end


%---
h.UserData = D;

%--- yoke other figures
if D.yoke && ~exist('foreignAction', 'var')
    h = findobj('TooltipString','coordinates (voxels)');
    h(h==D.hpos) = [];
    for i=1:length(h)
        h(i).String = D.hpos.String;
        h(i).Callback(h(i),[], true);
    end
end

end


