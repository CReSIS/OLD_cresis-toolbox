function [precisionType,nBytes] = TDMS_getTypeSeekSizes
%TDMS_getTypeSeekSizes  Returns information for fread and fseek
%
%   [precisionType,nBytes] = TDMS_getTypeSeekSizes
%
%   For internal use
%
%   precisionType : the string that goes to fread
%   nBytes        : gets multiplied by # of values to know how many bytes
%                   to skip when seeking
%

precisionType  = cell(1,10);
nBytes = zeros(1,68);

precisionType{1}  = 'int8=>int8';
nBytes(1) = 1;
precisionType{2}  = 'int16=>int16';
nBytes(2) = 2;
precisionType{3}  = 'int32=>int32';
nBytes(3) = 4;
precisionType{4}  = 'int64=>int64';
nBytes(4) = 8;
precisionType{5}  = 'uint8=>uint8';
nBytes(5) = 1;
precisionType{6}  = 'uint16=>uint16';
nBytes(6) = 2;
precisionType{7}  = 'uint32=>uint32';
nBytes(7) = 4;
precisionType{8}  = 'uint64=>uint64';
nBytes(8) = 8;
precisionType{9}  = 'single=>single';
nBytes(9) = 4;
precisionType{10} = 'double=>double';
nBytes(10) = 8;
nBytes(32) = 1; %Size is specified in bytes, not characters
nBytes(33) = 1;
nBytes(68) = 16;


end