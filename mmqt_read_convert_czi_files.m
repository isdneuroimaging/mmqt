function [fnamesCZI, fnamesMat, metadata] = mmqt_read_convert_czi_files(inCZI, summary, depth, orderChannels, figureScaling, figureHide)
% function [fnamesCZI, fnamesMat, metadata] = mmqt_read_convert_czi_files(dirCZI, summary, depth, orderChannels, figureScaling, figureHide)
% function [fnamesCZI, fnamesMat, metadata] = mmqt_read_convert_czi_files(fnameCZI, summary, depth, orderChannels, figureScaling, figureHide)
%
% Read Carl Zeiss Image (CZI) file, containing two color channels:
%   1. Layer: DAPI staining of nuclei
%   2. Layer: anti-Iba1 steining of microglia
%   (see also argument "orderChannels")
% In particularly, this script reads the the image data and the meta-data from the CZI-File.
% The image data will be converted into a 4D array. It saves the 4D matrix together with information
% about image size, image scaling and other meta-information into a Matlab (.mat) File, with the suffix "_stack1_raw.mat".
% The script will check if the found CZI-File have already been converted previously and skip them in this case.
%
% Input arguments:
%   fnameCZI        - path-name to single CZI-File, to be converted
%   dirCZI          - path-name to folder in which the CZI-File(s) are located
%   summary         - logical, indicating whether to show only the summary of the newly found CZI files. 
%                       If false, then in addition to the summary, newly found CZI files will be read and converted.
%                       (optional; defaults to summary=false)
%   depth           - number, indicating the search depth relative to the directory containing the CZI-files (i.e. dirCZI), as specified in file 'microglia_project_setup.m'
%                       If depth=0, all levels below dirCZI will be
%                       searchd.
%                       (optional; defaults to depth=1)
%   orderChannels   - vector of indices; order color channels according to this vector of ordinal numbers
%                       The order of channels should correspond to:
%                           1. DAPI
%                           2. anti-Iba
%                       If this is not the case, flip channels after reading from CZI-file, using this argument.
%                       Example: 
%                       If orderChannels=[2 1]; then the second color channel will become the first in the 4D matrix and vice versa
%                       If orderChannels=[], the order will not be changed
%                       IF orderChannels=[1 2], the order will not be changed, but if there was a third layer in the CZI, this layer will be removed
%                       (optional; defaults to orderChannels=[])
%   figureScaling   - 
%
% Output arguments:
%   fnamesCZI       - path-names of the found CZI-files; if input was a single CZI-file, then the input is returned unchanged
%   fnamesMat       - path-names of the created MAT-Files with image data (cell array of character vectors)
%   metadata        - meta-data read from the CZI-File(s), in table format

%% declare default values for unspecified or empty variables
varNames = {'summary', 'depth', 'orderChannels', 'figureScaling', 'figureHide'};
varDefault = {'false', '1', '[]', '2.5', 'false'};
for i=1:length(varNames)
    if ~exist(varNames{i}, 'var') || isempty(eval(varNames{i}))
        eval([varNames{i} '= ' varDefault{i} ';'])
    end
end

%% add Bio-Formats toolbox to the path (needed to read .czi files)
if ~exist('bfopen', 'file')
    %--- specify an absolute path to bfmatalb
    path_bfmatlab = '/Volumes/Local/Matlab/bfmatlab';
    if ~isdir(path_bfmatlab)
        %--- if not valid try to find bfmatalb in path relative to the running script
        dirname = fileparts(mfilename('fullpath'));
        path_bfmatlab = fullfile(dirname, 'bfmatlab');
        if ~isdir(path_bfmatlab)
            error('Bio-Formats toolbox not found! Make sure to add the toolbox to the Matlab search path!')
        end
    end
    addpath(path_bfmatlab);
end

%% add helper

%% check input
typeOfFile = exist(inCZI, 'file');
%--- regular file
if any(typeOfFile == [2 3 4 5 6])
    fnamesCZI = {inCZI};
    nCZI = 1;
elseif typeOfFile == 7
    %--- make sure there is no trailing slash in the dir path
    dirCZI = regexprep(inCZI , '\/$', '');
else
    error('Input is not a valid file or directory!')
end


%% search source folder for .czi files
%--- define search depth relative to dirCZI
if exist('dirCZI', 'var')
    if ~exist('depth','var') || ~isnumeric(depth)
        depth=1;
        fprintf('- NOTE: Search depth for CZI-files has been set to 1!\n')
    end
    if depth<=0
        depth = [filesep, '**', filesep];
    else
        depth = depth-1;
        depth = [filesep, repmat(['*' filesep],1,depth)];
    end
    %--- search CZI files
    fnamesCZI = glob([dirCZI, depth, '*.czi']);
    if isempty(fnamesCZI)
        error('Input directory does not contain any CZI File at the specified search depth')
    else
        %--- display file list
        fprintf('\nFound %g .czi files:\n', length(fnamesCZI))
        fprintf('  %s\n',fnamesCZI{:})
        nCZI = length(fnamesCZI);
        fprintf('(n=%d)\n\n', nCZI)
    end
end


%% check whether the CZI was read/converted already previously or was excluded (skip these files)
idxDone = false(length(fnamesCZI),1);
for iFn=1:length(fnamesCZI)
    %--- decompose path into dirname, basename and extension
    [folder, basename] = fileparts(fnamesCZI{iFn});
    %---
    fnameMatT = fullfile(folder, [basename '_stack1_raw.mat']);
    if exist(fnameMatT,'file')
        idxDone(iFn) = true;
    end
end
fnamesCZI(idxDone) = [];
if isempty(fnamesCZI)
    error('All files have been converted already previously')
end

%% Summarize (only if input was a directory)
if exist('dirCZI', 'var')
    fprintf('\nNewly detected CZI files:\n')
    fprintf('  %s\n', fnamesCZI{:})
    fprintf('\nSummary:\n')
    fprintf('N=%d files were found overall\n', nCZI)
    if exist('nExc','var')
        fprintf('N=%d files were excluded according to meta-data file\n', nExc)
    end
    if exist('nInc','var')
        fprintf('N=%d files were included according to meta-data file\n', nInc)
    end
    fprintf('N=%d files have been converted already previously\n', sum(idxDone))
    fprintf('N=%d files have been newly detected and will be converted\n', sum(~idxDone))
else
    fprintf('\nFollowing file has been provided as input and will be converted:\n  %s\n', fnamesCZI{1})
end
%% Return, if only the summary was requested
if exist('summary', 'var') && summary
    fnamesMat = [];
    metadata = [];
    return
end

%% loop through the list of czi files and import them into MAT-Files
ht = tic;
fnamesMat = cell(length(fnamesCZI),1);
for iFn=1:length(fnamesCZI)
        
    %% prepare folder/file-names for output
    %--- decompose path into dirname, basename and extension
    [folder, basename] = fileparts(fnamesCZI{iFn});
    %---
    fnamesMat{iFn} = fullfile(folder, [basename '_stack1_raw.mat']);
    fnamePng = fullfile(folder, [basename '_ortho1_raw.png']);
        
    %% Read the CZI file
    fprintf('\nReading %d. out of %d files:\n   %s\n', iFn, length(fnamesCZI), fnamesCZI{iFn})
    data = bfopen(fnamesCZI{iFn});
    %% Extract meta data, if possible
    try
        %% get metadata from CZI into a cell array, which we call "metaPairs"
        metaCZI = data{1, 2};
        metaCZI_iterator = metaCZI.keySet().iterator();
        metaPairs = cell(metaCZI.size(),2);
        for i=1:metaCZI.size()
            metaPairs{i,1} = metaCZI_iterator.nextElement();
            metaPairs{i,2} = metaCZI.get(metaPairs{i,1});
            % fprintf('%s = %s\n', pairs{i,:})
        end
        %% extract desired meta-data from CZI (extract wanted keys)
        keyWanted = {'Size[XYZC]','Scaling[XYZ]','Objective.*Model','PinholeDiameter','LaserPower'}; %,'Wavelength'};
        idx = [];
        metaPairsWanted = cell(0,2);
        for i=1:length(keyWanted)
            idx = find(~cellfun(@isempty,regexp(metaPairs(:,1), keyWanted{i})));
            temp = metaPairs(idx,:);
            for j=1:size(temp)
                if ~isnan(str2double(temp(j,2)))
                    temp{j,2} = str2double(temp(j,2));
                    if strcmp(keyWanted{i},'Scaling[XYZ]')
                        temp{j,2} = round(temp{j,2} *10^10) /10^4;
                    end
                end
            end
            metaPairsWanted = [metaPairsWanted; temp];
        end
        
        
        %% create a table "meta" with the desired meta-information
        varNames = regexprep(metaPairsWanted(:,1),'.*\|','');
        varNames = regexprep(varNames, ' *#1','');
        varNames = regexprep(varNames, ' *#','');
        % meta = cell2table([{fname, id}, pairsWanted(:,2)'],'VariableNames', [{'FileName'; 'ID'}; varNames(:)])
        meta = cell2table([{folder, basename}, metaPairsWanted(:,2)'],'VariableNames', [{'Folder'; 'Basename'}; varNames(:)]);
        
        %% add variables to the metadata-table
        meta.StackX = meta.SizeX * meta.ScalingX;
        meta.StackY = meta.SizeY * meta.ScalingY;
        meta.StackZ = meta.SizeZ * meta.ScalingZ;
        %---
        meta.include = 1;
        meta.dateConverted = {date};
        
    catch
        warning('%s\n         %s\n', 'Metadata could not be extracted!', 'Scaling of image will be assumed to be 1 micron in all directions!')
        meta = [cell2table({folder, basename},'VariableNames', {'Folder'; 'Basename'}),...
            array2table([nan(1, 4) ones(1,3)] , 'VariableNames', {'SizeX','SizeY','SizeZ','SizeC', 'ScalingX','ScalingY','ScalingZ'})];
    end

    %% Create image matrix from data structure
    %--- check meta data and, if incorrect, assume 2 color channels
    if size(data{1},1) ~= meta.SizeZ * meta.SizeC
        warning('%s\n         %s\n', 'Number of Z-planes and color-channels in the meta data does not correspond to size of data portion!',...
                'Assuming two color-channels and calculating the number of Z-planes accordingly!')
        meta.SizeC = 2;
        meta.SizeZ = size(data{1},1)/meta.SizeC;
        meta.SizeX = size(data{1}{1,1},1);
        meta.SizeY = size(data{1}{1,1},2);
    end
    %--- put data into 4D matrix
    k=0;
    img = zeros(meta{1,{'SizeX','SizeY','SizeZ','SizeC'}},'like',data{1}{1,1});
    for i=1:meta.SizeZ
        for j=1:meta.SizeC
            k=k+1;
            img(:,:,i,j) = data{1}{k,1};
        end
    end
    %--- order channels
    if ~isempty(orderChannels)
        if length(orderChannels) < meta.SizeC
            meta.SizeC = length(orderChannels);
        elseif length(orderChannels) > meta.SizeC
            orderChannels(meta.SizeC+1:end) = [];
        end
        img = img(:,:,:, orderChannels);
    end
    
    %% create header (with similar format as used for nifti files)
    clearvars hdr
    hdr.dim = [4 meta{1,{'SizeX','SizeY','SizeZ','SizeC'}}];
    hdr.pixdim = [4 meta{1,{'ScalingX','ScalingY','ScalingZ'}} 1];
    hdr.name = basename;
    fprintf('- Dimensions of image in micrometer:\n   %g %g %g\n', hdr.dim(2:4).*hdr.pixdim(2:4))
    %% show image in orthogonal view
    % close all
    h=show_cells_3D(img,hdr,[],[1 99],figureScaling, [], figureHide);
    
    %% save figure as PNG
    F = getframe(h);
    imwrite(F.cdata, fnamePng);    
    
    %% write image stack to mat-file
    fprintf('Saving stack to MAT-file:\n  %s\n', fnamesMat{iFn})
    save(fnamesMat{iFn}, 'img', 'hdr', 'meta')
    
    %% collect metadata for all images
    if iFn==1
        metadata = meta;
    else
        metadata = outerjoin(metadata, meta, 'MergeKeys', true);
    end
    %%
    display_time_delay
    
end

%%
if length(fnamesCZI)>1
    fprintf('\nFinished looping through images!\n')
    display_time_delay(ht)
end

