function h = plot_multiple_panels(varargin)
% function h = plot_multiple_panels(X, varLabel)
% function h = plot_multiple_panels(h, X, varLabel)
%
% Plot multiple variables in a tiled arangements of tiled axes
%
% Input:
%       h           - handle to figure; if missing, a new figure will be
%                     created (optional argument)
%       X           - matrix or table with columns corresponding to variables and
%                     rows to observations
%       varLabels   - cell array of strings, corresponding to the labels of
%                     the variables. Length of "varLabels" has to equal the
%                     number of columns in "X". If X is table, and varLabel
%                     is not provided, the names of table variables will be
%                     used instead.
%                     (optional argument)
%
% Output:
%       h           - handle to the created figure
%

%% get input arguments
varNames = {'h', 'X', 'varLabel'};
if ~ishandle(varargin{1})
    varNames(1) = [];
end
for i=1:length(varargin)
    eval([varNames{i} '= varargin{i};'])
end

%% number of variables
n = size(X,2);

%% check whether table or array, in case of table, variable names are used for titles, if input varLAbel is missing
if istable(X)
    if ~exist('varLabel','var') % || length(varLabel)~=n
        varLabel = X.Properties.VariableNames;
    end
    X = table2array(X);
    if size(X,2)>n
        error('Some table columns contain vectors instead of scalars!')
    end
end

%% create plots
if ~exist('h','var') || ~ishandle(h)
    h = figure;
end
nrows = ceil(sqrt(n));
ncols = ceil(n/nrows);
% xScatter = rand(size(X,1),1)/2+0.75;
xScatter = rand(size(X,1),1)/3+1-0.5/3;
for i=1:n
    ha = subplot(nrows, ncols, i);
    %--- boxplot
    boxplot(X(:,i))
    ha.XTickLabel = [];
    %--- add randomly scattered points
    hold on
    scatter(xScatter, X(:,i), 10, 'b', 'filled')
    %--- histogram
    % histogram(X(:,i),10)
    %---
    title(varLabel{i},'interpreter','none')
end
