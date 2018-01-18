function mmqt_select_cells_and_features(trunk, figureHide)
% function mmqt_select_cells_and_features(fnameCZI, figureHide)
% function mmqt_select_cells_and_features(fnameMat, figureHide)
%
% Function selects cells according to certain criteria (see below) and, in
% addition, it allows to keeps only shape features, according to a selection
% specified in the file "mmqt_feature_selection.xlsx" in column "include".
%
% Criteria for the selection of cells are:
% 1. Cells are touching the border of the Z-stack and the center of the
%   soma is closer than a certain threshold to this border. This threshold is
%   15 micrometer for the lateral borders of the stack (X/Y direction), and
%   8 micrometer for the bottom and top of the stack (Z direction).
% 2. Cells are touching each other and the centers of their soma are closer
%   together than 15 micrometer.
% 3. Cells whos soma contains 2 nuclei, which migth be due to the fact,
%   that the soma of different cells are touching each other.
%
% 
% Input arguments:
%   fnameCZI      - path to CZI-file; note however, that this path-name
%                   will only be used to reconstruct the path-name of the MAT-file with the raw data (see below).
%                   In order to be found, the correspinding MAT-file should be stored in the same folder as the CZI-file.
%   fnameMat      - path to MAT-file, which contains a 4D matrix with the raw image matrix and a
%                   structure with the information about image size and scaling.
%                   This MAT-file can be created by reading CZI-files using the script: mmqt_read_convert_czi_files.m
%   figureHide    - logical, indicating whether to make figures visible or not
%                   If figureHide=true, figures will be made invisible and used only for saving picturs of the screenshot;
%                   (optional, defaults to figureHide=fasle)


%%
fprintf('\n\nSummarizing results:\n\n')

%% decompose pathname
[folder, basename, ~] = fileparts(trunk);
basename = regexprep(basename, '_stack1_raw$', ''); %--- in case the MAT-file is given, remove the suffix, which is added by mmqt_read_convert_czi_files.m
trunk = fullfile(folder, basename);

%%
if ~exist('figureHide','var') || isempty(figureHide)
    figureHide = false;
end

%% read features
fnameFeatures = [trunk, '_features_of_cells.txt'];
D = readtable(fnameFeatures);


%% make a selection of good cells
%--- criterions are:
% 1. distance to the stack border 
% 2. distance to other interconnected cells
% 3. whether the soma resulted from merging several soma during skeletonization

%% definition of thresholds for distance criterions
thrDistBorderXY = 15;
thrDistBorderZ  = 8;
thrDistCells = 15;


%% extract distances of cells to the stack borders in all 6 directions and set distance=Inf, if a border is not touched
%--- index of columns containing the info, whether a cell touches the stack border in negative X/Y/Z and positive X/Y/Z direction
idxTouch = grep(fieldnames(D), 'stackTouch', 'l');
%--- index of columns containing the distance of the centroid of the soma to the the stack border in negative X/Y/Z and positive X/Y/Z direction
idxDist = grep(fieldnames(D), 'stackDistM', 'l');
%--- extract distances and substitute them with Inf, if the cell is not touching the border in a certain direction
distances = D{:,idxDist};
distances(D{:,idxTouch}==0) = Inf;

%% threshold the distance of cells to the border
distanceBorderOkXY = all(distances(:,[1 2 4 5]) > thrDistBorderXY, 2);
distanceBorderOkZ  = all(distances(:,[3 6]) > thrDistBorderZ, 2);
distanceBorderOk   = distanceBorderOkXY & distanceBorderOkZ;

%% threshold the distance between interconnected cells (i.e. distances of the centroid of the soma)
distanceCellsOk = D.groupShortestDist > thrDistCells | isnan(D.groupShortestDist);

%% logical index for good cells
D.cellOk = distanceBorderOk & distanceCellsOk & D.mergedSoma==0;
cellOkStr = repmat({'excluded'}, height(D), 1);
cellOkStr(D.cellOk) = {'included'};

fprintf('\nCells excluded/included:\n\n')
disp(crosstab_label(cellOkStr))


%% Adjust certain features
% %--- Adjust features which are calculated by deviding through the number of branches or branch segments:
% %--- if the number of branches is zero, the devision is not defined giving "NaN". These "NaN" values are here substituted by "0" to avoid loosing these cells for the analysis.
% varAdjust = grep(D.Properties.VariableNames,{'PerBranch', 'PerSegment'}, 't');
% if ~isempty(varAdjust)
%     D(D.branchN==0,varAdjust) = {0};
% end


%% read list of variables to be included
vars = readtable('mmqt_feature_selection.txt');
varNames = vars.variable(vars.include==1);

%--- Display names of included features
fprintf('Features/variables included:\n')
disp(varNames)

%--- include also the ID of the cell
varNames = unique(['soma'; varNames], 'stable');

%% keep only good cells and selected features in table
D = D(D.cellOk, varNames);
D.Properties.VariableNames{1} = 'ID_cell';


%% show summary
% summary(D(:,2:end))

%% plot features
h = figure;
h.Position = [45 5 1360, 1200];
if figureHide
    h.Visible = 'off';
end
if width(D) > 6
    plot_multiple_panels(h, D(:,2:end), D.Properties.VariableNames(2:end));
else
    [~, ax, ~] = gplotmatrix(D{:,2:end}, [], [],[],[],[], [], [], D.Properties.VariableNames(2:end));
    for i=1:numel(ax)
        ax(i).YLabel.Interpreter = 'none';
        ax(i).XLabel.Interpreter = 'none';
    end
end

%%
%--- save screenshot
F = getframe(h);
imwrite(F.cdata, [trunk,'_features_summarized.png'])


%% write the cell-wise statistics to Excel file
fnameFeaturesSelected = [trunk, '_features_summarized.txt'];
writetable(D, fnameFeaturesSelected, 'Delimiter', '\t');




