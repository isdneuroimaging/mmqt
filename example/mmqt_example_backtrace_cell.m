
%% Example shows how to backtrace a cell in the stacks using its ID
%
% In order to try this example, you need a dataset/stack alread being processed
% with MMQT. Change the pathname assigned to variable 'basename', to make
% it a valid pathname-trunk for your dataset/stack. The 'basename' will be 
% appended in this example with various suffix, corresponding to different 
% files created during processing with MMQT.
basename = '/Path/to/your/z-stack/z-stack-name-trunk';

%% Read the table with the features of each cell
% Please note that this table contains all cells, also those being
% excluded in the final table, summarizing cell features across all analysed stacks,
% if these cells are touching a border of the stack and if in addition the soma
% is located too close to a borders of the stack (cropped cells).
D = readtable([basename, '_features_of_cells.txt']);

%% Select a cell ID; for this example the cell with the minimum sphericity is selected
[minSphericit, id] = min(D.sphericity);

%% find the center of gravity of the cell soma
cog = round(D{id,5:7});

%% load the stacks
stackSeg = load([basename, '_stack5_segregated_cells.mat']);
stackEqu = load([basename, '_stack2_equalized.mat']);

%% crop the stack with equalized intensity histograms (the other stacks are already cropped)
crop = load([basename, '_prepro1_spatial_corrZ_sliceOk.mat']);
stackEqu.img = stackEqu.img(:,:,crop.idxSliceOk(1):crop.idxSliceOk(2),:);


%% visualize the stacks with equalized intensity histograms, with the segmentation, and with the segregation of cells
%--- the hair cross is positioned at the soma's center of gravity
show_cells_3D(stackSeg.img, stackSeg.hdr, cog)
show_cells_3D([basename, '_stack3_segmented.mat'],[], cog)
show_cells_3D(stackEqu.img, stackEqu.hdr, cog)

%% visualize only the selected cell
[x,y,z] = ind2sub(size(stackSeg.img), find(stackSeg.img==id));
imgCell = stackSeg.img(min(x):max(x), min(y):max(y), min(z):max(z));
imgCell(imgCell~=id)=0;
show_cells_3D(imgCell, stackSeg.hdr)

%% visualize maximum intensity projection
imgCellMax = max(imgCell,[],3);
figure;
imagesc(rot90(imgCellMax>0))
axis equal tight
colormap(gray)


