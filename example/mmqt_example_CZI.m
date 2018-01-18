function mmqt_example_CZI()
%
% This example shows how to process a single Carl Zeiss Image (CZI) file
% A valid CZI file, has to containing two color channels:
%   1. Layer: DAPI staining of nuclei
%   2. Layer: anti-Iba1 steining of microglia
% If the order of color layers is different in your own data, please see
% help of function mmqt_read_convert_czi_files.m for instructions



%% specify path to CZI file
%--- single file
fnameCZI = 'enter_here_the_path_to_the_example_CZI-file';
% e.g.: fnameCZI = '/Volumes/Local/microglia/data_czi/mouse_CB51/CB51_5dpdMCAO_40x_contra2.czi';

%% specify parameters for visualization
figureScaling = 2.5; %--- relevant only for the figures showing orthogonal slices through the stack; a value of one means that one pixel on the screen correponds to one pixel in the image
figureHide    = true; %--- creates figures invisibly, only for saving picturs of the screenshot;

%% Process the images

%--- start timer
ht = tic;

%--- read CZI file
mmqt_read_convert_czi_files(fnameCZI, [], [], [],figureScaling, figureHide);
%--- processing step 1
mmqt_segment_image(fnameCZI, figureScaling, figureHide);
%--- processing step 2
mmqt_skeletonize_microglia(fnameCZI);
%--- processing step 3
mmqt_extract_features(fnameCZI, figureScaling, figureHide);

%--- select cells and shape features
mmqt_select_cells_and_features(fnameCZI, figureHide)


fprintf('\n\nFinished looping through stacks!\n')
display_time_delay(ht)
