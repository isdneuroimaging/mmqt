function hdr = change_hdr_datatype_info(hdr,type)
% function hdr = change_hdr_datatype_info(hdr,type)
%
% Change data-type information in nifti header (hdr)!
% "type" can be given as a number or a string indicating the datatype.
% Possibilities are:
% number:  2        4        8        16         62
% string:  'uchar'  'int16'  'int32'  'float32'  'double'

%%
typeNum = [2,4,8,16,62];
typeStr = {'uchar','int16','int32','float32','double'};
bitpix = [8,16,32,32,64];

if ischar(type)
    idx = find(strcmpi(type,typeStr));
    if isempty(idx)
        error('Unsupported data type!\n');
    end
elseif isnumeric(type)
    idx = find(type==typeNum);
    if isempty(idx)
        error('Unsupported data type!\n');
    end
else
    error('Unsupported data type!\n');
end

hdr.datatype = typeNum(idx);
hdr.datatypestr = typeStr{idx};
hdr.bitpix = bitpix(idx);