function mmqt_example_MAT()
%
% This example shows how to process a single MATLAB formatted binary file 
% (MAT-file), containing the raw Z-stack data and meta information regarding 
% the scaling of the image. 
%
% PLEASE NOTE: This example is basically the same as that in mmqt_example_CZI.m
% It just leaves out the conversion from CZI to MAT format.
%
% If your own data are aquired in the Carl Zeiss Image (CZI) format, then
% you can possibly start directly from the CZI-files (see mmqt_example_CZI.m)
% Raw data in other formats have to be converted to a compatible MAT-file
% by the customer. The requirements for this MAT-file can be found in the
% README.txt file.%

%% specify path to CZI file
%--- single file
fnameMAT = 'enter_here_the_path_to_the_example_MAT-file';
% e.g.: fnameCZI = '/Volumes/Local/microglia/data_czi/mouse_CB51/CB51_5dpdMCAO_40x_contra2_stack1_raw.mat';

%% specify parameters for visualization
figureScaling = 2.5; %--- relevant only for the figures showing orthogonal slices through the stack; a value of one means that one pixel on the screen correponds to one pixel in the image
figureHide    = true; %--- creates figures invisibly, only for saving picturs of the screenshot;

%% Process the images

%--- start timer
ht = tic;



%--- processing step 1
mmqt_segment_image(fnameMAT, figureScaling, figureHide);
%--- processing step 2
mmqt_skeletonize_microglia(fnameMAT);
%--- processing step 3
mmqt_extract_features(fnameMAT, figureScaling, figureHide);

%--- select cells and shape features
mmqt_select_cells_and_features(fnameMAT, figureHide)


fprintf('\n\nFinished looping through stacks!\n')
display_time_delay(ht)
