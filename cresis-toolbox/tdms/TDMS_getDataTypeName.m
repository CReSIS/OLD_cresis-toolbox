function typeName = TDMS_getDataTypeName(dataType)
%TDMS_getDataTypeName  Returns a string indicating dataType
%
%   Given a Labview numeric datatype value, this returns a string
%
%   typeName = TDMS_getDataTypeName(dataType)

switch dataType
    case 0
        typeName = 'void';
    case 1
        typeName = 'int8';
    case 2
        typeName = 'int16';
    case 3
        typeName = 'int32';
    case 4
        typeName = 'int64';
    case 5
        typeName = 'uint8';
    case 6
        typeName = 'uint16';
    case 7
        typeName = 'uint32';
    case 8
        typeName = 'uint64';
    case 9
        typeName = 'single';
    case 10
        typeName = 'double';
    case 11
        typeName = 'extendedFloat';
    case 25
        typeName = 'singleWithUnit';
    case 26
        typeName = 'doubleWithUnit';
    case 27
        typeName = 'extendedFloatWithUnit';
    case 32
        typeName = 'string';
    case 33
        typeName = 'logical';
    case 68
        typeName = 'timestamp';
    case 2^32-1
         typeName = 'DAQmx';
    otherwise
        error('Unrecognized data type: %d',dataType)
end

end