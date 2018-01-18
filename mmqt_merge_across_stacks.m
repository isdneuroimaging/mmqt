function mmqt_merge_across_stacks(fnames, fnameOut)
% function mmqt_merge_across_stacks(fnames)
%
% Merges featrues, extracted from multiple Z-stacks
%

%%
fprintf('\n\nMerging features across stacks:\n\n')

for i=1:length(fnames)
    %% decompose pathname
    [folder, basename, ~] = fileparts(fnames{i});
    basename = regexprep(basename, '_stack1_raw$', ''); %--- in case the MAT-file is given, remove the suffix, which is added by mmqt_read_convert_czi_files.m
    trunk = fullfile(folder, basename);
    fnameFeatures = [trunk, '_features_summarized.txt'];
    
    
    %% read table
    D = readtable(fnameFeatures);
    %--- add column with name of stack (file basename)
    if exist('DD','var') && any(strcmp(DD.stack, basename))
        basename = [basename, '_', num2str(i)];
    end
    D.stack = repmat({basename},height(D),1);
    D = D(:,[end,1:end-1]);
    
    %% merge tables
    if i>1
        DD = outerjoin(DD, D, 'MergeKeys', true);
    else
        DD = D;
    end
    
end

%% group statistics for each stack
DDm = grpstats(DD(:,[1,3:end]), 'stack', 'median');

%%
if exist(fnameOut, 'file')
    error('Output file "%s" exists allready!\nPlease delete it first or choose different filename!', fnameOut)
end
[~,~,ext] = fileparts(fnameOut);
if any(strcmp(ext,{'.xlsx','.xls'}))
    writetable(DD, fnameOut, 'Sheet', 1);
    writetable(DDm, fnameOut, 'Sheet', 2);
    fprintf('Merged shape-feature tables were saved to: %s\n', fnameOut)
    fprintf('  Sheet1: features per cell\n')
    fprintf('  Sheet2: median of features per stack\n')
else
    writetable(DD, fnameOut);
    fprintf('Merged shape-feature tables were saved to: %s\n', fnameOut)
end

