function output = TDMS_initData(dataType,nSamples)
%TDMS_initData  Initializes raw data arrays
%
%   output = TDMS_initData(dataType,nSamples)
%
%   dataType - Labview datatype
%   nSamples - # of samples to initialize to
%
%   See Also: TDMS_getDataTypeName

switch dataType
    case 1
        output   = zeros(1,nSamples,'int8');
    case 2
        output   = zeros(1,nSamples,'int16');
    case 3
        output   = zeros(1,nSamples,'int32');
    case 4
        output   = zeros(1,nSamples,'int64');
    case 5
        output   = zeros(1,nSamples,'uint8');
    case 6
        output   = zeros(1,nSamples,'uint16');
    case 7
        output   = zeros(1,nSamples,'uint32');
    case 8
        output   = zeros(1,nSamples,'uint64');
    case {9 25}
        output   = zeros(1,nSamples,'single');
    case {10 26}
        output   = zeros(1,nSamples,'double');
    case 32
        output   = cell(1,nSamples); %string
    case 33
        output   = false(1,nSamples);
    case 68
        output   = zeros(1,nSamples,'double');
%     case intmax('uint32')
%         output   = zeros(1,nSamples,'int16');
    otherwise
        error('Unhandled data type for raw data: %s',TDMS_getDataTypeName(dataType))
end