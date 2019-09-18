function data = TDMS_readFileHelper_v1(fid,optionStruct,metaStruct,paramsStruct)
%TDMS_readFileHelper_v1
%
%

% Time stamps in TDMS are stored as a structure of two components:
% (i64) seconds: since the epoch 01/01/1904 00:00:00.00 UTC (using the Gregorian calendar and ignoring leap seconds)
% (u64) positive fractions: (2^-64) of a second
% Boolean values are stored as 1 byte each, where 1 represents TRUE and 0 represents FALSE.


SECONDS_IN_DAY  = paramsStruct.SECONDS_IN_DAY;
CONV_FACTOR     = paramsStruct.CONV_FACTOR;
UTC_DIFF        = paramsStruct.UTC_DIFF;
STRING_ENCODING = paramsStruct.STRING_ENCODING;

%INPUT UNPACKING
%==========================================================================
keepDataArray        = optionStruct.keepDataArray;
numValuesToGetActual = optionStruct.numValuesToGetActual;

rawDataInfo      = metaStruct.rawDataInfo;
segInfo          = metaStruct.segInfo;
numberDataPoints = metaStruct.numberDataPoints;

numObjects = length(rawDataInfo);
numSegs    = length(segInfo);

%INITIALIZATION OF OBJECTS
%==========================================================================
curFileIndex  = zeros(1,numObjects); %current # of samples read from file
curDataIndex  = zeros(1,numObjects); %current # of samples assigned to output
data          = cell(1,numObjects);  %a pointer for each channel
dataTypeArray = [rawDataInfo.dataType];

propNames  = cell(1,numObjects);
propValues = cell(1,numObjects);

for iObject = 1:numObjects
    propNames{iObject}  = rawDataInfo(iObject).propNames;
    propValues{iObject} = rawDataInfo(iObject).propValues;
    if numberDataPoints(iObject) > 0 && keepDataArray(iObject)
        data{iObject} = TDMS_initData(dataTypeArray(iObject),numberDataPoints(iObject));
    end
end

%==================================================================
%                        RAW DATA PROCESSSING
%==================================================================

%This will be used later for fread and fseek
%Simplifies the switch statements
[precisionType, nBytes] = TDMS_getTypeSeekSizes;

%Get end of file position, seek back to beginning
fseek(fid,0,1);
eofPosition = ftell(fid);
fseek(fid,0,-1);

for iSeg = 1:numSegs
    curSeg = segInfo(iSeg);
    %Seek to this raw position, this is needed to avoid meta data
    fseek(fid,curSeg.rawPos,'bof');
    
    nChunksUse = curSeg.nChunks;
    for iChunk = 1:nChunksUse
        %------------------------------------------------------------------
        %Interleaved data processing
        %------------------------------------------------------------------
        if curSeg.isInterleaved
            objOrder  = curSeg.objOrder;
            dataTypes = dataTypeArray(objOrder);
            nRead     = curSeg.nSamplesRead;
            
            %error checking
            if any(dataTypes ~= dataTypes(1))
                error('Interleaved data is assumed to be all of the same type')
            end
            
            if any(nRead ~= nRead(1))
                error('# of values to read are not all the same')
            end
            
            %NOTE: unlike below, these are arrays that we are working with
            startI = curFileIndex(objOrder) + 1;
            endI   = curFileIndex(objOrder) + nRead(1);
            curFileIndex(objOrder) = endI;
            curDataIndex(objOrder) = endI;
            
            nChans        = length(objOrder);
            numValuesRead = nRead(1);
            switch dataTypes(1)
                case {1 2 3 4 5 6 7 8 9 10}
                    temp = fread(fid,numValuesRead*nChans,precisionType{dataTypes(1)});
                case 32
                    error('Unexpected interleaved string data')
                    %In Labview 2009, the interleaved input is ignored
                    %Not sure about other versions
                case 33
                    %This never seems to be called, shows up as uint8 :/
                    temp = logical(fread(fid,numValuesRead*nChans,'uint8'));
                case 68
                    temp = fread(fid,numValuesRead*2*nChans,'*uint64');
                    temp = (double(temp(1:2:end))/2^64 + double(typecast(temp(2:2:end),'int64')))...
                        /SECONDS_IN_DAY + CONV_FACTOR + UTC_DIFF/24;
                otherwise
                    error('Unexpected data type: %d',dataTypes(1))
                    
            end
            
            %NOTE: When reshaping for interleaved, we must put nChans as
            %the rows, as that is the major indexing direction, we then
            %grab across columns
            %Channel 1 2 3 1  2  3
            %Data    1 2 3 11 22 33 becomes:
            %   1 11
            %   2 22    We can now grab rows to get individual channels
            %   3 33
            temp = reshape(temp,[nChans numValuesRead]);
            for iChan = 1:nChans
                if keepDataArray(objOrder(iChan))
                    data{objOrder(iChan)}(startI(iChan):endI(iChan)) = temp(iChan,:);
                end
            end
            
        else
            %--------------------------------------------------------------
            %NOT INTERLEAVED
            %--------------------------------------------------------------
            for iObjList = 1:length(curSeg.objOrder);
                I_object = curSeg.objOrder(iObjList);
                
                numValuesAvailable   = curSeg.nSamplesRead(iObjList);
                dataType             = dataTypeArray(I_object);
                
                curFileIndex(I_object) = curFileIndex(I_object) + numValuesAvailable;
                
                %Actual reading of data (or seeking past)
                %------------------------------------------
                if ~keepDataArray(I_object);
                    fseek(fid,numValuesAvailable*nBytes(dataType),'cof');
                else
                    startI = curDataIndex(I_object) + 1;
                    endI   = curDataIndex(I_object) + numValuesAvailable;
                    curDataIndex(I_object) = endI;
                    switch dataType
                        case {1 2 3 4 5 6 7 8 9 10}
                            data{I_object}(startI:endI) = fread(fid,numValuesAvailable,precisionType{dataType});
                        case 32
                            %Done above now ...
                            strOffsetArray = [0; fread(fid,numValuesAvailable,'uint32')];
                            offsetString = startI - 1;
                            for iString = 1:numValuesAvailable
                                offsetString = offsetString + 1;
                                temp = fread(fid,strOffsetArray(iString+1)-strOffsetArray(iString),'*uint8');
                                data{I_object}{offsetString}  = native2unicode(temp,STRING_ENCODING)'; %#ok<*N2UNI>
                            end
                            %NOTE: Even when using a subset, we
                            %will only ever have one valid read
                        case 33
                            data{I_object}(startI:endI)   = logical(fread(fid,numValuesAvailable,'uint8'));
                        case 68
                            temp = fread(fid,numValuesAvailable*2,'*uint64');
                            %First row: conversion to seconds
                            %Second row: conversion to days, and changing of reference frame
                            data{I_object}(startI:endI) = (double(temp(1:2:end))/2^64 + double(typecast(temp(2:2:end),'int64')))...
                                /SECONDS_IN_DAY + CONV_FACTOR + UTC_DIFF/24;
                        otherwise
                            error('Unexpected type: %d', dataType)
                    end
                end
            end
        end
    end
    
    
    %Some error checking just in case
    if iSeg ~= numSegs
        Ttag = fread(fid,1,'uint8');
        Dtag = fread(fid,1,'uint8');
        Stag = fread(fid,1,'uint8');
        mtag = fread(fid,1,'uint8');
        if ~(Ttag == 84 && Dtag == 68 && Stag == 83 && mtag == 109)
            error('Catastrophic error detected, code probably has an error somewhere')
        end
    else
        if eofPosition ~= ftell(fid) && ~metaStruct.eof_error
            error('Catastrophic error detected, code probably has an error somewhere')
        end
    end
end

%ERROR CHECKING ON # OF VALUES READ
%==========================================================================
if ~isequal(numValuesToGetActual,curDataIndex)
    error('The # of requested values does not equal the # of returned values, error in code likely')
end

%END OF READING RAW DATA
%==========================================================================
fclose(fid);
