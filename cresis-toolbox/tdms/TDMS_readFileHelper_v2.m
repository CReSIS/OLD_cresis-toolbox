function data = TDMS_readFileHelper_v2(fid,optionStruct,metaStruct,paramsStruct)
%TDMS_readFileHelper_v2
%
%

SECONDS_IN_DAY  = paramsStruct.SECONDS_IN_DAY;
CONV_FACTOR     = paramsStruct.CONV_FACTOR;
UTC_DIFF        = paramsStruct.UTC_DIFF;
STRING_ENCODING = paramsStruct.STRING_ENCODING;

%INPUT UNPACKING
%==========================================================================
subsetInfo           = optionStruct.subsetInfo;
numValuesToGetActual = optionStruct.numValuesToGetActual;

rawDataInfo      = metaStruct.rawDataInfo;
segInfo          = metaStruct.segInfo;

numObjects = length(rawDataInfo);

%INITIALIZATION OF OBJECTS
%==========================================================================
curDataIndex  = zeros(1,numObjects);
data          = cell(1,numObjects);  %a pointer for each channel
dataTypeArray = [rawDataInfo.dataType];

for iObject = 1:numObjects
    if numValuesToGetActual(iObject) > 0
        data{iObject} = TDMS_initData(dataTypeArray(iObject),numValuesToGetActual(iObject));
    end
end

%==================================================================
%                        RAW DATA PROCESSSING
%==================================================================
precisionType = TDMS_getTypeSeekSizes;

for iRead = 1:size(subsetInfo,1)
    
    %This code relies on subsetInfo, which gets defined in
    %TDMS_handleGetDataOption (originally called grabInfo)
    
    fseek(fid,subsetInfo(iRead,1),'bof');
    
    I_object      = subsetInfo(iRead,4);
    numValuesRead = subsetInfo(iRead,2);
    dataType      = dataTypeArray(I_object);
    
    startI = curDataIndex(I_object) + 1;
    endI   = curDataIndex(I_object) + numValuesRead;
    curDataIndex(I_object) = endI;
    switch dataType
        case {1 2 3 4 5 6 7 8 9 10} %numeric
            data{I_object}(startI:endI) = fread(fid,numValuesRead,precisionType{dataType});
        case 32 %strings
            curSeg             = segInfo(subsetInfo(iRead,5));
            numValuesAvailable = curSeg.nSamplesRead(curSeg.objOrder == I_object);
            strOffsetArray     = [0; fread(fid,numValuesAvailable,'uint32')];
            if subsetInfo(iRead,3) ~= 1
               %see TDMS_handleGetDataOption for more info on this
               %need to do some seeking
               fseek(fid,strOffsetArray(subsetInfo(iRead,3)),'cof');
            end
            offsetString = startI - 1;
            for iString = subsetInfo(iRead,3):subsetInfo(iRead,3)+subsetInfo(iRead,2)-1
                offsetString = offsetString + 1;
                temp         = fread(fid,strOffsetArray(iString+1)-strOffsetArray(iString),'*uint8');
                data{I_object}{offsetString} = native2unicode(temp,STRING_ENCODING)'; %#ok<*N2UNI>
            end
        case 33 %boolean
            data{I_object}(startI:endI) = logical(fread(fid,numValuesRead,'uint8'));
        case 68 %timestamps
            temp = fread(fid,numValuesRead*2,'*uint64');
            %First row: conversion to seconds
            %Second row: conversion to days, and changing of reference frame
            data{I_object}(startI:endI) = (double(temp(1:2:end))/2^64 + double(typecast(temp(2:2:end),'int64')))...
                /SECONDS_IN_DAY + CONV_FACTOR + UTC_DIFF/24;
        otherwise
            error('Unexpected type: %d', dataType)
    end
end

%ERROR CHECKING ON # OF VALUES READ
%==========================================================================
if ~isequal(numValuesToGetActual,curDataIndex)
    error('The # of requested values does not equal the # of returned values, error in code likely')
end
%numValuesToGetActual vs curDataIndex

%END OF READING RAW DATA
%==========================================================================
fclose(fid);
