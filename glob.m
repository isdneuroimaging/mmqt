function varargout = glob(pattern, type, parse)
% function [pathnames, date] = glob(pattern, type, parse)
%
% Input arguments:
%   pattern - Character vector or cell array of character vectors, specifying each one search pattern.
%             This argument is basically the same as the "name" argument to the "dir" function (not a regular expression)
%             However, in contrast to the "dir" function, multiple patterns can be specified as a cell array of character vectors (i.e. strings) 
%   type    - Character specifying the file type to look for:
%               - f regular file
%               - d directory
%             by default both, directories and regular files, will be returned
%   parse   - Integer indicating whether to parse pathname. Depending on the
%             value, either the whole path, the dirname or the basename will be returned.
%               - If argument is missing or empty, the whole path will be returned.
%               - If parse==1 the dirname
%               - If parse==2 the basename
%               - If parse==3 the basename without extension
%               - If parse==4 the dirnam basename extension in separate colums of a cell array
%
% Output arguments:
%   pathnames - Cell array of character strings containing the full path-names of the found files
%   date      - Cell array of character strings containing the modification date
%
% example patterns:
% pattern = '*/*'; %--- all files/folders in subdirectory of depth 2
% pattern = '*/*/*'; %--- all files/folders in subdirectory of depth 3
% pattern = '**/*'; %--- recursive search of all files/folders relative to the current folder
% pattern = '/Volumes/projects/**/*'; %--- recursive search of all files/folders relative to the folder /Volumes/projects/
% pattern = '*analyse*.m' %--- all files with extension .m and containing the word analyse, in the current directory
% pattern = '**/*analyse*.m' %--- all files with extension .m and containing the word analyse, searching recursively, starting in the current directory
%


%% call the dir command to search the file hierarchy
if ~iscell(pattern) && ischar(pattern)
    L = dir(pattern);
else
    L = struct([]);
    for i=1:length(pattern)
        L = [L; dir(pattern{i})];
    end
end

%%
%% remove "." and ".." file placeholder
L(ismember({L.name}, {'.','..','.DS_Store'})) = [];

%% select file type
if exist('type','var') && ~isempty(type) && ischar(type)
    switch type
        case {'d', '-d'}
            L = L([L.isdir], :);
        case {'f', '-f'}
            L = L(~[L.isdir], :);
    end
end
    

%%
for i=1:length(L)
   [~,L(i).nameX,L(i).ext] = fileparts(L(i).name);
end

%% create output arguments
if isempty(L) %--- avoid problems with non-existing fields (i.e. nameX, ext)
    varargout{1} = {};
elseif exist('parse','var') && ~isempty(parse) && isnumeric(parse)
    switch parse
        case 1
            varargout{1} = {L.folder}';
        case 2
            varargout{1} = {L.name}';
        case 3
            varargout{1} = {L.nameX}';
        case 4
            varargout{1} = [{L.folder}',{L.nameX}',{L.ext}'];
        otherwise
            %--- create full pathname (by concatenating dirname and basename)
            varargout{1} = fullfile({L.folder}, {L.name})';
    end
else
    %--- create full pathname (by concatenating dirname and basename)
    varargout{1} = fullfile({L.folder}, {L.name})';
end


if nargout > 1
    varargout{2} = {L.date}';
end




