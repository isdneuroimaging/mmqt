function hdr = change_hdr_dim(hdr,sizeImg)
% function hdr = change_hdr_dim(hdr,sizeImg)

%--- remove ones from the end of the vector
while ~isempty(sizeImg) && sizeImg(end)==1
    sizeImg(end) = [];
end

if length(sizeImg)<3
    hdr.dim(1:4) = [3 1 1 1];
    hdr.dim(5:end) = 0;
    hdr.dim(2:length(sizeImg)+1) = sizeImg;
    hdr.pixdim(5:end) = 0;
else
    hdr.dim(1:length(sizeImg)+1) = [length(sizeImg) sizeImg];
    hdr.dim(length(sizeImg)+2:end) = 0;
    hdr.pixdim(length(sizeImg)+2:end) = 0;
end

