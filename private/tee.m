function varargout = tee(fid, formatSpec, varargin)
% function varargout = tee(fileID, formatSpec, A1,...,An)
% function varargout = tee(fileName, formatSpec, A1,...,An)
% function varargout = tee([], formatSpec, A1,...,An)
%
% This function is an extension to fprintf/sprintf, which allows to duplicate 
% the formated string in one of three ways:
% 1. File & Command Window
% 2. File & Output Variable
% 3. Output Variable & Command Window
%
% Details: 
% With a valid file ID or file-name as first input and without output argument, 
% this function writes the formated text to the file & the Command Window.
% If an output argument is provided, then the formated text is NOT written
% to the command window, but returned in the output variable instead.
% If the first input argument is empty, the output is always written to the
% Command Window and optionally an output argument (if provided).
%
% First input argument has to be a file ID or a file name (path name).
% If an file name is given, then the file is opened with permissions 'a+'
% and closed after the writing operation.
% All further input arguments are expected to be as in fprintf

%% vector of logicals defining actions to perform
actions = true(1,3); %--- [file, varargout, CommandWindow]

%% verify validity of first input argument
if isempty(fid)
    actions(1) = false;
    closeFid = false;
elseif isnumeric(fid) && any(fid == fopen('all'))
    closeFid = false;
elseif ischar(fid)
    if ~isempty(regexp(fid, '[,*%?"<>|]', 'once'))
        error('First input has to be a valid file identifier or a valid file name!')
    else
        fid = fopen(fid, 'a+');
        closeFid = true;
    end
else
    error('First input has to be a valid file identifier, a valid file name, or empty (i.e. [])!')
end

%% check whether an output argument is desired
if nargout==0
    actions(2) = false;
end

%% decide whether to print to Command Window
if actions(1) && actions(2)
    actions(3) = false;
end

%% write the formated string
%--- write to file
if actions(1)
    fprintf(fid, formatSpec, varargin{:});
end
%--- write to output argument
if actions(2)
    varargout{1} = sprintf(formatSpec, varargin{:});
    for k = 2:nargout
        varargout{k} = [];
    end
end

%--- write to Command Window
if actions(3)
    fprintf(formatSpec, varargin{:});
end

%% close the file, if it was opened at the beginning
if closeFid
    fclose(fid);
end

