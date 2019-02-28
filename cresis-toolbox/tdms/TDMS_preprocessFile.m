function metaStruct = TDMS_preprocessFile(fid,tdmsFileName,params)
%TDMS_preprocessFile  Populates props and gets segment information
%
%   NOTE: This file shouldn't be called directly by the user
%
%   metaStruct = TDMS_preprocessFile(fid,tdmsFileName,params)
%
%   INPUTS
%   =======================================================================
%   fid          : file id of opened tdms file, may be null if INDEX_DEBUG is true
%   tdmsFileName : full path to tdms file OR tdms index file for INDEX_DEBUG
%   params       : all optional inputs from TDMS_readTDMSFile
%
%   OUTPUTS
%   =======================================================================
%   metaStruct :
%       structure with fields:
%
%   The following 3 variables are all the same length, and have one index
%   for each unique object.  Each object may be one of the following types,
%   file root, group object, or channel object.  Only channel objects may
%   have data.
%
%   rawDataInfo           : (structure array)
%           NOTE: This holds the most recent info for an object, and can be
%           changed on a per segment basis, segInfo is used to hold
%           information that is needed per segement
%
%           .lengthOfIndex    - length of raw data index (not used)
%           .dataType         - enumerated dataType (Labview)
%           .dimensionData    - currently ignored (should be 1d)
%           .numberOfValues   - # of values to read each time
%           .totalSizeBytes   - # of bytes per values
%           .numberProperties - # of properties defined for the object
%           .infoSet          - a logical on whether or not this data has
%                               been set (for error checking)
%           .propNames        - a cell array for definining property names
%           .propValues       - cell array for holding property values
%
%   numberDataPoints     : an array which specifies the # of data points for
%                          each object (used for preallocation)
%   objectNameList       : (cell array of strings) names of all objects in the
%                          TDMS file
%
%   Info for each segment ...
%
%   segInfo              : (struct array)
%           .rawPos - position of raw data for that segment
%           .kTocNewObjList - flag on whether a new channel list was
%           created or not
%           .objOrder - array specifying indices of the objects to assign
%                       data to, order specifies read order, and values
%                       are indexed
%           into the objectNameList & rawDataInfo struct
%           .nRawObjects - length of objOrder
%           .nSamplesRead - array, # of samples to read for each object
%           .isInterleaved - logical, whether that segment uses interleaved
%                               data or not
%           .isBigEndian - logical
%           .nChunks     - # of times to read data using the parameters
%                           above before moving onto the next segment
%
%   See also: TDMS_processLeadIn, TDMS_getDataSize, TDMS_getPropValue


%PARAMETERS - passed in from TDMS_readTDMSFile
%==========================================================================
UTC_DIFF        = params.UTC_DIFF;
MAX_NUM_OBJECTS = params.MAX_NUM_OBJECTS;
MAX_NUM_PROPS   = params.MAX_NUM_PROPS;
N_SEGS_GUESS    = params.N_SEGS_GUESS;
N_SEGS_INC      = params.N_SEGS_INC;
DEBUG           = params.DEBUG;
DATE_STR_FORMAT = params.DATE_STR_FORMAT;
UNICODE_FORMAT  = params.STRING_ENCODING;
USE_INDEX       = params.USE_INDEX;
INIT_CHUNK_SIZE = params.INIT_CHUNK_SIZE;
TDMS_INDEX_EXT  = params.TDMS_INDEX_EXT;
MACHINE_FORMAT  = params.MACHINE_FORMAT;
STRING_ENCODING = params.STRING_ENCODING;
CURRENT_VERSION = params.CURRENT_VERSION;
INDEX_DEBUG     = params.INDEX_DEBUG;

%DONT CHANGE THESE
%==========================================================================
LEAD_IN_LENGTH = 28; %# of bytes of lead in

%FIGURING OUT WHICH FILE TO READ
%==========================================================================
[tdmsPathToFile,tdmsNameOnly] = fileparts(tdmsFileName);
if isempty(params.META_STRUCT)
    isIndexFID = false;
    if USE_INDEX || INDEX_DEBUG
        if INDEX_DEBUG
            indexFile = tdmsFileName;
        else
            indexFile = fullfile(tdmsPathToFile,[tdmsNameOnly TDMS_INDEX_EXT]);
        end
        if exist(indexFile,'file')
            fid = fopen(indexFile,'r',MACHINE_FORMAT,STRING_ENCODING);
            isIndexFID = true;
        else
            %Just use the tdms file, which we have already tested to exist
            if INDEX_DEBUG
                %NOTE: With INDEX_DEBUG we explicitly passed in the index
                %file to read and parse (generally for debugging purposes
                %so if it doesn't exist, then we have a problem
                error('Specified tdms_index file doesn not exist')
            end
        end
    end
else
    %Quits early if meta struct input is valid
    %======================================================================
    metaStruct = params.META_STRUCT;
    if ~isstruct(metaStruct) || ~isfield(metaStruct,'version') || ~isfield(metaStruct,'fileName')
        error('The META_STRUCT parameter passed should be a structure with fields "version" and "fileName"')
    end
    if metaStruct.version ~= CURRENT_VERSION
        error('Mismatch in version creation, version of struct %d, please run with latest version: %d')
    end
    if ~strcmpi(metaStruct.fileName,tdmsNameOnly)
        fprintf('Filename from metaStruct: %s\n',metaStruct.fileName);
        fprintf('Filename from input: %s\n',tdmsNameOnly);
        error('Mismatch in filenames, see printout above')
    end
    return
end

if isIndexFID
    lastLetter = double('h');  %used for .tdms_index files
else
    lastLetter = double('m');  %used for .tdms files
end


%INITIALIZE OUTPUTS
%==========================================================================
curNumTotalObjects  = 0; %Current # of objects
numberDataPoints    = zeros(1,MAX_NUM_OBJECTS);
objectNameList      = cell(1,MAX_NUM_OBJECTS);  %names of all objects

rawDataInfo = initRawInfoStruct(MAX_NUM_PROPS,MAX_NUM_OBJECTS);
segInfo     = initSegStruct(N_SEGS_GUESS);
%==========================================================================

%TEMP VARIABLES
%==========================================================================
objectHasRawData    = false(1,MAX_NUM_OBJECTS); %This is needed for adding
%# of data points together
ranOnce    = false;
nSegs      = 0;

%Get eof & return to start
fseek(fid,0,1);
eofPosition = ftell(fid);
fseek(fid,0,-1);

%START OF READING
%==========================================================================
curPosCounter = 0; %Used to keep track of where in the actual tdms file we 
%would be at various points, this is used for parsing verification as well
%as assisting with more complicated reads
while ftell(fid) ~= eofPosition
    
    nSegs = nSegs + 1;
    if nSegs > length(segInfo)
        segInfo = [segInfo initSegStruct(N_SEGS_INC)]; %#ok<AGROW>
    end
    
    if DEBUG
        disp('------  START OF DEBUG -------')
        fprintf(2,'CURRENT SEGMENT: %d\n',nSegs);
        fprintf(2,'Current file position: %d\n',curPosCounter);
    end
    
    %LEAD IN HANDLING
    %======================================================================
    [flags,info,eof_error]  = TDMS_processLeadIn(fid,lastLetter);
    
    if eof_error
        %This should only happen once at the end 
        fprintf(2,['WARNING: File was not closed properly.\n' ... 
            'Data will most likely be missing at the end of the file\n']);
        nSegs = nSegs - 1;
        break
    end
    
    curPosCounter = curPosCounter + LEAD_IN_LENGTH + info.segLength;
    segInfo(nSegs).rawPos = curPosCounter - info.segLength + info.metaLength;
    
    
    if ~ranOnce
        %NOTE: This might be false if no channels are defined ...
        if flags.kTocNewObjList == false
            error('the kTocNewObjList was false when first run, this is not expected')
        end
        ranOnce = true;
    end
    %======================================================================
    
    
    %======================================================================
    %                       META INFORMATION PROCESSING
    %======================================================================
    if flags.hasMetaData
        
        %Get # of changed objects
        numNewObjInSeg = fread(fid,1,'uint32');
        
        %Reinitialize order list if new
        if flags.kTocNewObjList
            curObjListCount      = 0;
            objOrder             = zeros(1,2*numNewObjInSeg);
            nValuesRead          = zeros(1,2*numNewObjInSeg);
            %NOTE: I padded this by doubling the #, we might append extra
            %channels in subsequent reads, I currently don't resize this
            %...
        end
        
        for iNewObject = 1:numNewObjInSeg
            
            %1) GET OBJECT PATH
            %-------------------------------
            objPathLength = fread(fid,1,'uint32');
            temp          = fread(fid,objPathLength,'*uint8');
            objName       = native2unicode(temp,UNICODE_FORMAT)'; %#ok<*N2UNI>
            
            %POPULATE ALL OBJECT LIST
            %----------------------------------------
            objIndex = find(strcmp(objectNameList(1:curNumTotalObjects),objName),1);
            if isempty(objIndex)
                curNumTotalObjects  = curNumTotalObjects + 1;
                if curNumTotalObjects > length(rawDataInfo)
                    rawDataInfo      = [rawDataInfo      ...
                        initRawInfoStruct(MAX_NUM_PROPS,MAX_NUM_OBJECTS)]; %#ok<AGROW>
                    numberDataPoints = [numberDataPoints zeros(1,MAX_NUM_OBJECTS)]; %#ok<AGROW>
                    objectNameList   = [objectNameList   cell(1,MAX_NUM_OBJECTS)]; %#ok<AGROW>
                    objectHasRawData = [objectHasRawData false(1,MAX_NUM_OBJECTS)];  %#ok<AGROW>
                end
                objIndex            = curNumTotalObjects;
                objectNameList{curNumTotalObjects} = objName;
            end
            
            
            
            %RAW DATA INDEX PROCESSING
            %---------------------------------------------
            curPos = ftell(fid);
            
            rawDataIndexLength = fread(fid,1,'uint32');
            
            rawDataInfo(objIndex).lengthOfIndex = rawDataIndexLength;
            
            if DEBUG
                fprintf(2,'RawDataLength: %d\n',rawDataIndexLength);
                fprintf(2,'CurrentPos: %d\n',curPos);
            end
            
            
            switch rawDataIndexLength
                case 0 %Same as previous
                    if rawDataInfo(objIndex).infoSet == false
                        error('Channel %s set to use previous rawDataIndex but this channel is new',objName)
                    end
                    
                    %NOTE: "same as previous segment" apparently means
                    %"same as the previous one with data ..."
                    
                    if rawDataInfo(objIndex).numberOfValues > 0 && ~objectHasRawData(objIndex)
                        objectHasRawData(objIndex) = true;
                    end
                    
                case 2^32-1 %no raw data
                    rawDataInfo(objIndex).infoSet    = true;
                    objectHasRawData(objIndex)       = false;
                otherwise
                    objectHasRawData(objIndex)    = true;
                    rawDataInfo(objIndex).infoSet = true;
                    
                    %DATA TYPE HANDLING
                    %------------------------------------------------------
                    dataType = fread(fid,1,'uint32');
                    if rawDataInfo(objIndex).infoSet && dataType ~= rawDataInfo(objIndex).dataType && numberDataPoints(objIndex) > 0
                        error('Raw data type for channel %s has changed from %d to %d',...
                            objName,rawDataInfo(objIndex).dataType,dataType)
                    else
                        rawDataInfo(objIndex).dataType = dataType;
                    end

                    %DATA SIZE HANDLING
                    %-----------------------------------------------------
                    rawDataInfo(objIndex).dimensionData  = fread(fid,1,'uint32');
                    if rawDataInfo(objIndex).dimensionData ~= 1
                        error('Code doesn''t yet handle non 1D data')
                    end
                    
                    rawDataInfo(objIndex).numberOfValues = fread(fid,1,'uint64');
                    
                    %RawDaqMX
                    %------------------------------------------------------
                    %NOTES:
                    if dataType == 2^32-1
                        %JIM CODE IN PROGRESS
                        rawDataInfo(objIndex).isRawDAQmx = true; %We can
% % %                         %post process this to convert from bits to an
% % %                         %actual value
% % %                         
% % %                         %FORMAT:
% % %                         %1) -
% % %                         daqMXVersion = rawDataIndexLength;
% % %                         if ~ismember(daqMXVersion,[4713 4714])
% % %                             error('Unexpected version: %d',daqMXVersion)
% % %                         end
% % % 
% % % 
% % %                         %2) Let's get the remaining # of bytes
% % %                         %This seems to be 32 ...
% % %                         
% % % % %                         wtf = fread(fid,8,'uint32');
% % % % %                         disp(wtf)
% % %                         %4713
% % %                         %13000012
% % %                         %15000018
% % %                         %? -> 8th indicates data type? OR -> bytes per sample
% % %                         
% % %                         %4714
% % %                         %1 0 0 0 256 256 512
% % % 
% % %                         %3) Update datetype - where the heck is it?
% % %                         %rawDataInfo(objIndex).dataType = dataType;
                    else %Non RawDAQmx type
                        if rawDataInfo(objIndex).dataType == 32
                            %If string, size is specified by an additional field
                            rawDataInfo(objIndex).totalSizeBytes = fread(fid,1,'uint64');
                        else
                            rawDataInfo(objIndex).totalSizeBytes = ...
                                rawDataInfo(objIndex).numberOfValues*TDMS_getDataSize(dataType);
                        end
                        
                        %Another chance to check correct reading
                        if curPos + rawDataIndexLength ~= ftell(fid)
                            error(['Raw Data Index length was incorrect: %d stated vs %d observed,' ...
                                'likely indicates bad code or a bad tdms file'],rawDataIndexLength,ftell(fid) - curPos)
                        end
                    end
                    
                    if DEBUG
                        fprintf(2,'nSegs: %d\n',nSegs);
                        fprintf(2,'objName: %s\n',objName);
                    end
            end
            
            
            %--------------------------------------------------------------
            %POPULATE ORDER TO RETRIEVE RAW DATA
            %--------------------------------------------------------------
            if objectHasRawData(objIndex)
                appendToList = false;
                if flags.kTocNewObjList
                    appendToList = true;
                else %Only append if not currently specified
                    I_objOrder = find(objOrder(1:curObjListCount) == objIndex,1);
                    if isempty(I_objOrder)
                        appendToList = true;
                    else
                        nValuesRead(I_objOrder) = rawDataInfo(objIndex).numberOfValues;
                    end
                end
                
                %NOTE: No overflow code in place yet, however we do
                %initialize with twice the # of objects specified to be 
                %in a segement when a new list is created, new segments
                %might add more objects to the list
                if appendToList
                    curObjListCount                 = curObjListCount + 1;
                    objOrder(curObjListCount)       = objIndex;
                    nValuesRead(curObjListCount)    = rawDataInfo(objIndex).numberOfValues;
                end
            end
            
            
            %PROPERTY HANDLING
            %--------------------------------------------------------------
            numberProperties = fread(fid,1,'uint32');
            %Below is the # of props already assigned to that channel
            nPropsChan       = rawDataInfo(objIndex).numberProperties;
            curProps         = rawDataInfo(objIndex).propNames;
            
            for iProp = 1:numberProperties
                propNameLength  = fread(fid,1,'uint32');
                temp            = fread(fid,propNameLength,'*uint8');
                propName        = native2unicode(temp,UNICODE_FORMAT)';
                propDataType    = fread(fid,1,'uint32');
                
                propIndex = find(strcmp(curProps(1:nPropsChan),propName),1);
                if isempty(propIndex)
                    %Updates needed for new properties
                    nPropsChan              = nPropsChan + 1;
                    propIndex               = nPropsChan;
                    curProps{propIndex}     = propName;
                    rawDataInfo(objIndex).propNames{propIndex} = propName;
                end
               
                %Update value
                propValue = TDMS_getPropValue(fid,propDataType,UTC_DIFF,DATE_STR_FORMAT);
                rawDataInfo(objIndex).propValues{propIndex} = propValue;
            end
            rawDataInfo(objIndex).numberProperties = nPropsChan;
            
            if DEBUG
                fprintf(2,'end of index position: %d\n',ftell(fid));
            end
        end
    end
    %======================================================================
    %               END OF META DATA PROCESSING
    
    
    
    %RAW DATA SAMPLE COUNTING
    %======================================================================
    
    segInfo(nSegs).kTocNewObjList = flags.kTocNewObjList;
    byteSizeRaw = info.segLength - info.metaLength;
    if ~flags.hasRawData || byteSizeRaw == 0
        segInfo(nSegs).nChunks = 0;
    else
        segInfo(nSegs).objOrder      = objOrder(1:curObjListCount);
        segInfo(nSegs).nRawObjects   = curObjListCount;
        segInfo(nSegs).nSamplesRead  = nValuesRead(1:curObjListCount);
        segInfo(nSegs).isInterleaved = flags.isInterleaved;
        segInfo(nSegs).isBigEndian   = flags.isBigEndian;
        
        %# OF CHUNK PROCESSING
        %------------------------------------------------------------------
        totalBytesPerChunk = sum([rawDataInfo(objOrder(1:curObjListCount)).totalSizeBytes]);
        nChunks            = byteSizeRaw/totalBytesPerChunk;
        
        %Some error checking
        %------------------------------------------
        if DEBUG
            fprintf(2,'nChunks: %d\n',nChunks);
            fprintf(2,'nSamplesRead: %s\n',mat2str(nValuesRead(1:curObjListCount)));
            fprintf(2,'totalBytesPerChunk: %d\n',totalBytesPerChunk);
            fprintf(2,'byteSizeRaw: %d\n',byteSizeRaw);
        end
        
        if nChunks ~= floor(nChunks)
            error(['The remaining data doesn''t split evently into' ...
                ' chunks, estimated # of chunks: %d'],nChunks)
        end
        
        chunkByteOffset = 0;
        %Increment the number of data points
        for iObject = 1:curObjListCount
            
            curIndex = objOrder(iObject);
            
            nSamplesReadCurObject = rawDataInfo(curIndex).numberOfValues;
            
            %Apparently some writers don't stripe the objects when some
            %objects don't have any data and some do ...
            if nSamplesReadCurObject > 0            
                %This allows us to grow these values if we haven't sufficiently preallocated
                if rawDataInfo(curIndex).chunkIndex + nChunks > rawDataInfo(curIndex).chunkLength
                    rawDataInfo(curIndex).dataMatrix  = [rawDataInfo(curIndex).dataMatrix; zeros(INIT_CHUNK_SIZE,3)];
                    rawDataInfo(curIndex).chunkLength = rawDataInfo(curIndex).chunkLength + INIT_CHUNK_SIZE;
                end

                %DataMatrix:
                %==============================================================
                %This information is used for reading parts of an object during
                %a single read instead of the entire object

                %Column 1, file position
                indices = (rawDataInfo(curIndex).chunkIndex+1):(rawDataInfo(curIndex).chunkIndex+nChunks);
                rawDataInfo(curIndex).dataMatrix(indices,1) = ...
                    segInfo(nSegs).rawPos + chunkByteOffset + (0:totalBytesPerChunk:(nChunks-1)*totalBytesPerChunk);

                %Column 2, first sample at that position
                rawDataInfo(curIndex).dataMatrix(indices,2) = ...
                    numberDataPoints(curIndex) + (0:nSamplesReadCurObject:(nChunks-1)*nSamplesReadCurObject) + 1;

                %Column 3, segment number
                rawDataInfo(curIndex).dataMatrix(indices,3) = nSegs;

                chunkByteOffset                  = chunkByteOffset + rawDataInfo(curIndex).totalSizeBytes;
                numberDataPoints(curIndex)       = numberDataPoints(curIndex) + nSamplesReadCurObject*nChunks;
                rawDataInfo(curIndex).chunkIndex = rawDataInfo(curIndex).chunkIndex+nChunks;
            end
        end
        
        %nChunksAll = nChunksAll + nChunks;
        segInfo(nSegs).nChunks = nChunks;
        
        %This needs to be handled in some manner for RawDaqMx
        if ~isIndexFID
            fseek(fid,byteSizeRaw,'cof');
        end
    end
    
    %Addition for raw data (needs to be fixed)
    %In general, the tdms_index requires skipping some meta data
    %I think this is only needed for raw data (I need to flush this out
    %better as it will cause an error when reading index files for regular
    %tdms files)
%     if isIndexFID
%         fseek(fid,segInfo(nSegs).rawPos,'bof');
%     end
end





%Trim output:
%==========================================================================
rawDataInfo         = rawDataInfo(1:curNumTotalObjects);
for iObject = 1:length(rawDataInfo)
    nProps = rawDataInfo(iObject).numberProperties;
    rawDataInfo(iObject).propNames  = rawDataInfo(iObject).propNames(1:nProps);
    rawDataInfo(iObject).propValues = rawDataInfo(iObject).propValues(1:nProps);
    
    %Add extra value for later processing
    if rawDataInfo(iObject).chunkIndex ~= 0
        rawDataInfo(iObject).chunkIndex = rawDataInfo(iObject).chunkIndex + 1;
        rawDataInfo(iObject).dataMatrix(rawDataInfo(iObject).chunkIndex,2) = numberDataPoints(iObject) + 1;
    end
end

numberDataPoints    = numberDataPoints(1:curNumTotalObjects);
objectNameList     = objectNameList(1:curNumTotalObjects);
segInfo             = segInfo(1:nSegs);

metaStruct = struct(...
    'eof_error',        eof_error,...
    'numberDataPoints', numberDataPoints,...
    'objectNameList',   {objectNameList},...
    'segInfo',          segInfo,...
    'rawDataInfo',      rawDataInfo,...
    'fileName',         tdmsNameOnly,...
    'version',          CURRENT_VERSION);

if isIndexFID
    %If using the index file to parse meta, close it
    fclose(fid);
end


end

function segStruct = initSegStruct(nSegs)
%initSegStruct
%
%   segStruct = initSegStruct(nSegs)
%
segStruct = struct(...
    'rawPos',repmat({0},[1 nSegs]),...
    'kTocNewObjList',0,...
    'objOrder',[],...
    'nRawObjects',0,...
    'nSamplesRead',[],...
    'isInterleaved',false,...
    'isBigEndian',false,...
    'nChunks',false);

end

function rawInfoStruct = initRawInfoStruct(MAX_NUM_PROPS,MAX_NUM_OBJECTS)

rawInfoStruct = struct( ...
    'isRawDAQmx',       false, ...
    'lengthOfIndex',    0,...
    'dataType',         0,...
    'dimensionData',    0,...
    'numberOfValues',   0,...
    'totalSizeBytes',   0,... %Only valid for strings
    'numberProperties', 0,...
    'chunkIndex',       0,...
    ... %Could initialize to zero, this would save space
    ... %with non-raw data objects
    'chunkLength',      0,...
    ... %SEE ALSO: chunk loop above for resizing
    'dataMatrix',       zeros(0,3),...
    ...%Position #
    ...%1st Sample #
    ...%Seg #
    'infoSet',          false, ...
    'propNames',        repmat({cell(1,MAX_NUM_PROPS)},1,MAX_NUM_OBJECTS),...
    'propValues',       repmat({cell(1,MAX_NUM_PROPS)},1,MAX_NUM_OBJECTS));
end
