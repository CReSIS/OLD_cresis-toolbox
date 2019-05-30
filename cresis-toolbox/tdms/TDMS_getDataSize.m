function dataSize = TDMS_getDataSize(dataType)
%TDMS_getDataSize  Returns data size in bytes
%
%   This function is used to return the size of raw data for predicting the
%   # of chunks of raw data
%
%   dataSize = TDMS_getDataSize(dataType)
%
%   INPUTS
%   =================
%   dataType : Labview dataType
%   
%   OUTPUTS
%   ==================
%   dataSize : size in bytes 
%
%   CALLED BY:
%   - TDMS_preprocessFile
%   - TDMS_handleGetDataOption


switch dataType
    case 1 %int8
        dataSize = 1;
    case 2 %int16
        dataSize = 2;
    case 3 %int32
        dataSize = 4;
    case 4 %int64
        dataSize = 8;
    case 5 %uint8
        dataSize = 1;
    case 6 %uint16
        dataSize = 2;
    case 7 %uint32
        dataSize = 4;
    case 8 %uint64
        dataSize = 8;
    case 9 %Single
        dataSize = 4;
    case 10 %Double
        dataSize = 8;
    case 25 %Single with unit
        dataSize = 4;
    case 26 %Double with unit
        dataSize = 8;
    case 32
        error('The size of strings is variable, this shouldn''t be called')
    case 33 %logical
        dataSize = 1;
    case 68 %timestamp => uint64, int64
        dataSize = 16;
%     case intmax('uint32')
%         %DAQmx
%         dataSize = 2; %Will need to be changed
%         %keyboard
    otherwise
        switch dataType
            case 0
                unhandledType = 'Void';
            case 11
                unhandledType = 'Extended float';
                %SIZE: 12 bytes
            case 27
                unhandledType = 'Extended float with unit';
                %SIZE: 12 bytes
            case intmax('uint32')
                unhandledType = 'DAQmx';
                %SIZE: Sent email about this
        end
        error('Unhandled property type: %s',unhandledType)
        %IMPROVEMENT:
        %We could fail silently and document this in the
        %structure (how to read DAQmx (how big to skip?)
end