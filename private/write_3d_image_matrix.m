function write_3d_image_matrix(hdr, img, fname)
% function write_3d_image_matrix(hdr, img, fname)
%
% This function is currently used only to allow an easy replacement of the
% previously used function write_nifti_gz
% In the future this function might be modified to allow saving in SRC
% compatible format which can be read by the DSI studio
%%
save(fname, 'hdr', 'img')