function propValue = TDMS_getPropValue(fid,propDataType,UTC_DIFF,DATE_STR_FORMAT)
%TDMS_getPropValue  Returns the property value given the Labview DataType
%
%   propValue = TDMS_getPropValue(fid,propDataType,UTC_DIFF)
%
%   See Also: TDMS_getDataTypeName

SECONDS_IN_DAY  = 86400;
CONV_FACTOR     = 695422; %datenum('01-Jan-1904')
UNICODE_FORMAT  = 'UTF-8';

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
switch propDataType
    case 1
        propValue    = fread(fid,1,'*int8');
    case 2
        propValue    = fread(fid,1,'*int16');
    case 3
        propValue    = fread(fid,1,'*int32');
    case 4
        propValue    = fread(fid,1,'*int64');
    case 5
        propValue    = fread(fid,1,'*uint8');
    case 6
        propValue    = fread(fid,1,'*uint16');
    case 7
        propValue    = fread(fid,1,'*uint32');
    case 8
        propValue    = fread(fid,1,'*uint64');
    case 9
        propValue    = fread(fid,1,'*single');
    case 10
        propValue    = fread(fid,1,'*double');
    case 25
        propValue    = fread(fid,1,'*single');
    case 26
        propValue    = fread(fid,1,'*double');
    case 32
        stringLength = fread(fid,1,'uint32');
        temp         = fread(fid,stringLength,'*uint8');
        propValue    = native2unicode(temp,UNICODE_FORMAT)';  %#ok<*N2UNI>
        %propValue    = fread(fid,stringLength,'*char')';
    case 33
        propValue = logical(fread(fid,1,'*uint8'));
    case 68
        %Eeek, uint64 really????
        %Time in seconds since 01/01/1904 0 UTC
        %Matlab: days since 01/00/0000 0 i.e. 1 represents 01/01/0000
        firstByte  = fread(fid,1,'uint64');
        secondByte = fread(fid,1,'int64');
        tSeconds     = firstByte/(2^64)+secondByte;
        propValue    = datestr(tSeconds/SECONDS_IN_DAY + CONV_FACTOR + UTC_DIFF/24,DATE_STR_FORMAT);
    otherwise
        error('Unhandled property type: %s',TDMS_getDataTypeName(propDataType))
end