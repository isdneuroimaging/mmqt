function [hdr, img, other] = read_3D_image_matrix(fname)
% function [hdr, img, other] = read_3D_image_matrix(fname)
%
% This function loads a Matlab MAT-file and tries to identify the variables
% containing the image matrix and the header.
% In the future this function might be modified to allow reading SRC
% files which can be read by the DSI studio

%% load the mat-file
S = load(fname);

%% find the field with the image matrix and the header
%--- initialize variables
img  = [];
hdr  = [];
other= [];
%--- go through fields of the loaded structure
ff = fieldnames(S);
for i=1:length(ff)
    if isnumeric(S.(ff{i})) && (length(size(S.(ff{i})))==3 || length(size(S.(ff{i})))==4)
        img = cast(S.(ff{i}), 'double'); %--- NOTE: to be compatible with read_nifti, the data will always be output as double
    elseif isstruct(S.(ff{i})) && isfield(S.(ff{i}),'pixdim')
        hdr = S.(ff{i});
    else
        other.(ff{i}) = S.(ff{i}); 
    end
end
