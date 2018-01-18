function varargout = grep(names, pattern, options)
% function lines = grep(filenames, pattern [, options])
% function lines = grep(lines, pattern [, options])
%
% Function works similar as the bash functions "grep" and "cat", depending on the chosen options: 
%
% It finds all cells in the cell array of character-strings given in 
% "filenames" or "lines", containing a certain pattern. 
% Alternatively, it can be used to first read the text in the files specified 
% in "filenames", and then it finds all lines in the read text, containing 
% a certain "pattern".
% The "pattern" argument can contain one or multiple patterns to search for
% in the input, and will return then all cells/lines, containing at least one of these patterns.
%
% The function can also be used to check the existence of files.
%
% For all options see description of "options" below.
%
% Input:
%       filenames - character vector or cell array with names of files containing text.
%       lines     - character vector or cell array with character strings.
%       pattern   - character vector or cell array with regular expression(s) to be looked for in each line of the input
%       options   - (optional argument) option flag character(s) specifying one ore more options:
%                       Option flags indicating what to do, if "filenames" are provided as input:
%                       'f'   read files and apply grep to the content (the default action for filenames)
%                       't'   treat filenames as text and apply grep directly on the filenames
%                       Further implemented options:
%                       'i'   ignore case in the regular expression (default is case-sensitve)
%                       'n'   negate the result; return the lines of the text, not containing the regular expression in "pattern"!
%                       'l'   return vector of type logical, indexing the text-line containing the regular expression in "pattern"
%                       'x'   Find files in "filenames" containing the "pattern", and return the names of these files or (with option 'l') a vector of logicals indexing these files. For multiple occurrences, each file is listed once only
%                             This option entails automatically option 'f' to be used.
%                             If option 'l' is selected a, vector of logicals with length equal length(filenames) will be returned.
%                       'e'   Test if files in "filenames" exist; returns names of existing files or (with option 'l') a vector of logicals indexing these files.
%                             This option entails automatically option 'f' to be used. Options 't', 'x' and 'i' will be ignored.
%                             The input argument 'pattern' will be ignored, but has to be given.
%                       'd'   As 'e' but tests only for directories
%                       'r'   As 'e' but tests only for regular files
%                       's'   Silence warning messages specific to this function (subfunctions might still throw warning messages)
%
% Output:
%       lines - cell array containing all cells in "filenames"/"lines" (or 
%               all lines read from the file(s) in "filenames"), in which the pattern was found
%
%
% Example:
%       %--- example returning vector of logials, being true for each cell in testNames, not containing the testPattern
%       testNames = {'34msyhHalloerg', 'asdgw40HalloWrt', 'wer 547Halla', 'Drtwef_f5Hallo', 'Drtwef7Halla'};
%       testPattern = '[^0-9][0-9]Hall[oa]$';
%       lines = grep(testNames, testPattern, 'nli') 
%       %--- example testing whether files in the "filenames" arguemnt exist
%       testNames = {'/Volumes/Local/projects/test.txt', '/Volumes/Local/projects/documentation.xlsx'}
%       grep(testNames, [], 'el')

%% check options
if ~exist('options','var') || ~ischar(options)
    options = '';
end
action = regexprep(options, '[^ft]', '');
ignoreCase = ~isempty(regexprep(options, '[^i]', ''));
negateResult = ~isempty(regexprep(options, '[^n]', ''));
returnLogical = ~isempty(regexprep(options, '[^l]', ''));
returnFilenames = ~isempty(regexprep(options, '[^x]', ''));
if ~isempty(regexp(options,'r', 'once'))
    checkExistence = [2 3 4 5 6];
elseif ~isempty(regexp(options,'d', 'once'))
    checkExistence = 7;
elseif ~isempty(regexp(options,'e', 'once'))
    checkExistence = [2 3 4 5 6 7];
else
    checkExistence = [];
end
% checkExistence = ~isempty(regexprep(options, '[^e]', ''));
waringOff = ~isempty(regexprep(options, '[^s]', ''));

%%
if waringOff
    warning('off', 'MATLAB:GREP:COSTUM')
else
    warning('on', 'MATLAB:GREP:COSTUM')
end
%% make sure the first 2 input arguments come as cell arrays of strings
if isempty(names)
    warning('MATLAB:GREP:COSTUM','First input argument was empty!')
    if returnLogical
        
        varargout{1} = false;
    else
        varargout{1} = [];
    end
    return
end
if ~iscell(names) && ischar(names)
    names = {names};
end
if isempty(pattern) && ~ischar(pattern) && isempty(checkExistence)
    error('Second input argument (regular-expression) not specified!')
%     warning('MATLAB:GREP:COSTUM','Second input argument (regular-expression) not specified!')
%     if returnLogical
%         varargout{1} = false(size(names));
%     else
%         varargout{1} = [];
%     end
end
if ~iscell(pattern) && ischar(pattern)
    pattern = {pattern};
end

%% spacial treatment of option 'e'
if ~isempty(checkExistence)
    idxFound = false(length(names),1);
    for i=1:length(names)
        if ismember(exist(names{i}, 'file'), checkExistence)
            idxFound(i) = true;
        end
    end
    if negateResult
        idxFound = ~idxFound;
    end
    if ~returnLogical
        idxFound = names(idxFound);
    end 
    if nargout>0
        varargout{1} = idxFound;
        for k = 2:nargout
            varargout{k} = [];
        end
    else
        if isempty(idxFound)
            fprintf('\n')
            fprintf('the resulting set of file-names is empty!\n')
        else
            if returnLogical
                temp = {'false','true'};
                idxFound = temp(idxFound+1);
            end
            fprintf('\n')
            fprintf(' %s\n',idxFound{:});
        end
    end
    return
end


%% check whether the first input contains file-names or character strings
readFiles = true;
if ~returnFilenames
    if isempty(action)
        %--- check whether all strings in names are valid file names
        for i=1:length(names)
            if ismember(exist(names{i}, 'file'),[0, 7]) %--- file does not exist or is a folder
                readFiles = false;
                break
            end
        end
    elseif strcmp(action,'t')
        readFiles = false;
    end
end
%% read lines from files
if readFiles
    kk = zeros(length(names),1);
    linesCat = cell(0);
    k=0;
    for i=1:length(names)
        if ~ismember(exist(names{i}, 'file'),[0, 7])
            fid = fopen(names{i});
            tline = fgetl(fid);
            while ischar(tline)
                k=k+1;
                linesCat{k,1} = tline;
                tline = fgetl(fid);
            end
            
            fclose(fid);
        else
            warning('MATLAB:GREP:COSTUM','File does not exist: %s', names{i})
        end
        kk(i) = k;
    end
    %---
    if k==0 %--- no text could be read from files
        fprintf('\n')
        warning('MATLAB:GREP:COSTUM','files were all either non-existing or empty!')
        if ~returnLogical
            if nargout>0
                varargout{1} = linesCat;
                for k = 2:nargout
                    varargout{k} = [];
                end
            end
            return
        end
    end
else
    linesCat = names;
end

%% apply grep to lines of text
%--- for testing
% linesCat = {'34msyhHalloerg', 'asdgw40HalloWrt', 'wer 547Halla', 'Drtwef_f5Hallo', 'Drtwef7Halla'};
% pattern = '[^0-9][0-9]Hall[oa]$';
%--- define anonymous function with placeholder for the pattern
if ignoreCase
    afunStr = 'afun = @(x) ~isempty(regexpi(x, ''%s''));';
else
    afunStr = 'afun = @(x) ~isempty(regexp(x, ''%s''));';
end
%--- apply anonymous function for each pattern
idxFound = false(size(linesCat));
patternFound = cell(size(linesCat));
idxPattern = nan(size(linesCat));
for i=1:length(pattern)
    if strcmp(pattern{i}, '') %--- allow to use an empty string as a pattern
        idxFoundT = strcmp(linesCat, '');
    else
        eval(sprintf(afunStr, pattern{i}));
        idxFoundT = cellfun(afun, linesCat);
    end
    if any(idxFoundT)
        % patternFound{idxFoundT} = pattern{i};
        patternFound(idxFoundT) = pattern(i);
        idxPattern(idxFoundT) = i;
    end
    idxFound = idxFound | idxFoundT;
end



%% consider further options
%---
if negateResult && ~returnFilenames
    idxFound = ~idxFound;
end


%---
if returnFilenames
    idxFoundT = find(idxFound);
    kk = [0; kk];
    for i=1:length(kk)-1
        if any(idxFoundT > kk(i) & idxFoundT <= kk(i+1))
            kk(i) = 1;
        else
            kk(i) = 0;
        end
    end
    idxFound = kk(1:end-1)>0;
    if negateResult
        idxFound = ~idxFound;
    end
end

%---
if returnLogical
    linesCat = idxFound;
else
    if returnFilenames
        linesCat = names(idxFound);
    else
        linesCat = linesCat(idxFound);
        patternFound = patternFound(idxFound);
        idxPattern = idxPattern(idxFound);
    end
end


%% if output requested, return the lines as output, otherwise display in command window
if nargout>0
    varargout{1} = linesCat;
    kk=2;
    if nargout>=2
        varargout{2} = patternFound;
        kk=3;
    end
    if nargout==3
        varargout{3} = idxPattern;
        kk=4;
    end
    for k = kk:nargout
        varargout{k} = [];
    end
else
    if isempty(linesCat)
        fprintf('\n')
        warning('MATLAB:GREP:COSTUM','the resulting set of lines is empty!')
    else
        if returnLogical
            temp = {'false','true'};
            linesCat = temp(linesCat+1);
        end
        fprintf('\n')
        fprintf(' %s\n',linesCat{:});
    end
end
