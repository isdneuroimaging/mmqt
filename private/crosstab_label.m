function [M, chi2, p] = crosstab_label(varargin)
% [M, chi2, p] = crosstab_label(Var1, Var2, options)
%
% Create Cross-tabulation with function "crosstab" and in addition add
% labels to the rows and columns of the table.
% Input:
%       Var1    - First variable
%       Var2    - Second variable (optional argument) 
%       options - Option flag character(s) specifying one ore more options (optional argument):
%                   's'   add summs along rows and colums to the table
%                   (only one option implemented yet)
%                              
% Output:
%       M       - Table with cross-tabulation and labels
%       chi2    - same as output from function "crosstab"
%       p       - same as output from function "crosstab"
%
%%
if ischar(varargin{end})
    options = varargin{end};
    varargin = varargin(1:end-1);
    nn = length(varargin);
else
    options = '';
    nn = nargin;
end
%---
addSums = regexprep(options, '[^s]', '');
%%
if nn==1
    [tbl, chi2, p, labels] = crosstab(varargin{1});
    M = array2table(tbl,'RowNames',labels(:,1),'VariableNames',{'Counts'});
elseif nn==2
    [tbl, chi2, p, labels] = crosstab(varargin{1},varargin{2});
    LT = ~cellfun(@isempty, labels);
    %%
    varNames = labels(LT(:,2),2);
    varNames = regexprep(varNames,'[^a-zA-Z0-9_]','_');
    varNames = regexprep(varNames,'^[0-9_]','v$0');
    %%
    M = array2table(tbl,'RowNames',labels(LT(:,1),1),'VariableNames', varNames);
else
    error('number of input arguments has to be 1 or 2')
end

%%
if ~isempty(addSums)
    M(end+1, :) = num2cell(sum(M{:,:},1));
    M(:, end+1) = num2cell(sum(M{:,:},2));
    M.Properties.RowNames{end} = 'sums_columns';
    M.Properties.VariableNames{end} = 'sums_rows';
end
