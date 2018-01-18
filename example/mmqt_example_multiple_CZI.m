function mmqt_example_multiple_CZI()
%
% This example shows how multiple CZI-files could be processed using a
% loop. Please notice, that the example data-set contains only one CZI-file.
% Therefore, this example will just go through one iteration of the for loop.
% However, if you want to make sure it really works with multiple
% CZI-files, you could just duplicate the example CZI-file for this test.


%% specify main folder, containing subfolder with CZI files (note: only one subfolder/CZI-file provided in example dataset)

%--- specify main folder here
directory_main = 'enter_here_the_path_to_the_main_data_folder';
% e.g.: directory_main = '/Volumes/Local/microglia/data_czi/'; %--- this folder should contain the folder "/mouse_CB51"

%--- find the CZI-files, included in subfolder of the main folder, using glob
fnamesCZI = glob([directory_main, '/*/*.czi']);
fprintf('Found following CZI-files:\n')
disp(fnamesCZI)

%% specify parameters for visualization
figureScaling = 2.5; %--- relevant only for the figures showing orthogonal slices through the stack; a value of one means that one pixel on the screen corresponds to one pixel in the image
figureHide    = true; %--- creates figures invisibly, only for saving picturs of the screenshot;

%% Process the images
%--- start timer the for loop
ht = tic;

for i=1:length(fnamesCZI)
    
    fprintf('\nProcessing %d. out of %d stacks\n', i, length(fnamesCZI))
    
    %--- read CZI file
    mmqt_read_convert_czi_files(fnamesCZI{i}, [], [], [],figureScaling, figureHide);
    
    %--- processing step 1
    mmqt_segment_image(fnamesCZI{i}, figureScaling, figureHide);
    %--- processing step 2
    mmqt_skeletonize_microglia(fnamesCZI{i});
    %--- processing step 3
    mmqt_extract_features(fnamesCZI{i}, figureScaling, figureHide);
    
    %--- select cells and shape features
    mmqt_select_cells_and_features(fnamesCZI{i}, figureHide)

    %--- close figures (which might be invisible) and diplay time needed for this loop
    close all
    fprintf('\nFinished %d. loop\n', i)
    display_time_delay(ht)
    
end

%--- time duration for the whole loop
fprintf('\n\nFinished looping through stacks!\n')
display_time_delay(ht)


%% merging features across CZI-Files
fnameOut = 'merged_features.xlsx';
%- Please note: if Excel is not availabel, save as speardsheet
%- (e.g. fnameOut = 'merged_features.csv')
mmqt_merge_across_stacks(fnamesCZI, fnameOut)

