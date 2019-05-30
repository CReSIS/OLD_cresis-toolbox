function varargout = TDMS_handleGetDataOption(mode,optionsIn,metaStruct)
%TDMS_handleGetDataOption  Processes input options 
%
%   Shouldn't be called directly by user
%
%   This function does some initial verfication on the types of inputs
%   that have been passed in.  The second call generates instructions
%   on how to read the file, given the meta data.
%
%   INITIAL CALL:
%   TDMS_handleGetDataOption('check',optionsIn)
%
%   SECOND CALL:
%   optionStruct = TDMS_handleGetDataOption('getArray',optionsIn,metaStruct)
%
%   See Also:
%       TDMS_retrievingSubsets

GET_DATA_OPTION  = optionsIn.GET_DATA_OPTION;
OBJECTS_GET      = optionsIn.OBJECTS_GET;
OBJECTS_IGNORE   = optionsIn.OBJECTS_IGNORE;
SUBSET_GET       = optionsIn.SUBSET_GET;
GET_INDICES      = optionsIn.GET_INDICES;
SUBSET_IS_LENGTH = optionsIn.SUBSET_IS_LENGTH;

if strcmp(mode,'check')
    
    %KEEP DATA ARRAY HANDLING
    %================================================
    switch lower(GET_DATA_OPTION)
        case {'getall' 'getnone'}
            if ~isempty(OBJECTS_GET) || ~isempty(OBJECTS_IGNORE)
                error('For GET_DATA_OPTION: %s, neither OBJECTS_GET or OBJECTS_IGNORE should be specified',GET_DATA_OPTION)
            end
        case 'getsubset'
            if isempty(OBJECTS_GET) || ~isempty(OBJECTS_IGNORE)
                error('For GET_DATA_OPTION: %s, OBJECTS_GET should be specified, OBJECTS_IGNORE should not',GET_DATA_OPTION)
            end
            if length(OBJECTS_GET) ~= 1
                error('For GET_DATA_OPTION: %s, OBJECTS_GET should have length 1, might be missing {} in struct initialization',GET_DATA_OPTION)
            end
            if ~isfield(OBJECTS_GET,'fullPathsKeep') && ~isfield(OBJECTS_GET,'groupsKeep')
                error('OBJECTS_GET must have fields fullPathsKeep and/or groupsKeep')
            end
            
        case 'ignoresubset'
            if ~isempty(OBJECTS_GET) || isempty(OBJECTS_IGNORE)
                error('For GET_DATA_OPTION: %s, OBJECTS_GET should be specified, OBJECTS_IGNORE should not',GET_DATA_OPTION)
            end
            if length(OBJECTS_IGNORE) ~= 1
                error('For GET_DATA_OPTION: %s, OBJECTS_IGNORE should have length 1, might be missing {} in struct initialization',GET_DATA_OPTION)
            end
            if ~isfield(OBJECTS_IGNORE,'fullPathsIgnore') && ~isfield(OBJECTS_IGNORE,'groupsIgnore')
                error('OBJECTS_GET must have fields fullPathsIgnore and/or groupsIgnore')
            end
        otherwise
            error('Unrecognized GET_DATA_OPTION')
    end
    
    %SUBSET RULE HANDLING
    %================================================
    if ~isempty(SUBSET_GET)
        if length(SUBSET_GET) ~= 2
            error(['SUBSET_GET is an array with columns' ...
                ' [startIndex Length], # columns observed: %d'],size(SUBSET_GET,2))
        end
        if SUBSET_GET(1) < 1
            error('The first input to SUBSET_GET should be 1 or higher')
        end
        if SUBSET_GET(2) < 1
            error('The second input to SUBSET_GET should be 1 or higher')
        end
        if strcmpi('getnone',GET_DATA_OPTION)
            error('The SUBSET_GET input option shouldn''t be used in conjunction with GET_DATA_OPTION of getNone')
        end
        
        %Could add additional check on 2 vs 1 if SUBSET_IS_LENGTH is false
    end
    
    %GET INDICES
    %==================================================
    if ~isempty(GET_INDICES)
        if ~strcmpi('getall',GET_DATA_OPTION)
            error('The GET_DATA_OPTION isn''t used with GET_INDICES, and should be left as getAll')
        end 
        if ~isstruct(GET_INDICES) || ~isfield(GET_INDICES,'group') || ...
                ~isfield(GET_INDICES,'channel') || ~isfield(GET_INDICES,'indices')
           error(['GET_INDICES must be a structure that defines the following fields: ' ...
               'group, channel, & indices'])
        end
        if ~isempty(SUBSET_GET)
            error('SUBSET_GET shouldn''t be defined when GET_INDICES is used')
        end
        %skipping for now, index checking ...
    end
    
else
    
    %metaStruct unpacking
    %-----------------------------------------
    objectPaths = metaStruct.objectNameList; 
    groupNames  = metaStruct.groupNames;
    rawDataInfo = metaStruct.rawDataInfo; %Information on a per object basis
    segInfo     = metaStruct.segInfo; %Information on a per segment basis
    
    numObjects = length(objectPaths);
    %KEEP DATA ARRAY HANDLING
    %======================================================================
    %NOTE: This variable is ultimately ignored if GET_INDICES is specified
    switch lower(GET_DATA_OPTION)
        case 'getall';
            keepDataArray = true(1,numObjects);
        case 'getsubset'
            if ~isfield(OBJECTS_GET,'fullPathsKeep')
                OBJECTS_GET.fullPathsKeep = {};
            end
            if ~isfield(OBJECTS_GET,'groupsKeep')
                OBJECTS_GET.groupsKeep = {};
            end
            keepDataArray = ismember(objectPaths,OBJECTS_GET.fullPathsKeep) | ...
                ismember(groupNames,OBJECTS_GET.groupsKeep);
        case 'ignoresubset'
            if ~isfield(OBJECTS_IGNORE,'fullPathsIgnore')
                OBJECTS_IGNORE.fullPathsIgnore = {};
            end
            if ~isfield(OBJECTS_IGNORE,'groupsIgnore')
                OBJECTS_IGNORE.groupsIgnore = {};
            end
            keepDataArray = ~(ismember(objectPaths,OBJECTS_IGNORE.fullPathsIgnore) | ...
                ismember(groupNames,OBJECTS_IGNORE.groupsIgnore));
        case 'getnone'
            keepDataArray = false(1,numObjects);
    end
    %Result, keepDataArray specifies whether or not to keep the data
    %for each object
   
    %SUBSET HANDLING
    %======================================================================
    if isempty(SUBSET_GET) && isempty(GET_INDICES)
        useSubset  = false;
        subsetInfo = struct([]);
        numValuesToGetActual = metaStruct.numberDataPoints;
        numValuesToGetActual(~keepDataArray) = 0;
    else
        useSubset = true;
        
        %In this section we create the indexStruct structure array.
        %The Length of indexStruct is equal to the # of raw data objects
        %we'll retrieve raw data from
        %------------------------------------------------------------------
        %.id      - which object in the file the entry refers to
        %.indices - [startIndex stopIndex]  OR [startIndex GrabLength]
        %   indices may have multiple rows, corresponding to multiple subsets
        if ~isempty(SUBSET_GET)
            %Here we keep objects based on keepDataArray
            objectsUse = find((metaStruct.numberDataPoints > 0 & keepDataArray));
            indexStruct = struct('id',num2cell(objectsUse),'indices',[SUBSET_GET(1) SUBSET_GET(2)]);
        else
            %-> path format -> /'group'/'channel'
            %this variable with contain the full paths of the objects
            %to retrieve
            objPaths_getIndex = arrayfun(@(x) ['/''' x.group '''/''' x.channel ''''],GET_INDICES,'UniformOutput',false);
            [~,loc] = ismember(objPaths_getIndex,objectPaths);
            if any(loc == 0)
                disp('Bad Paths:')
                disp(objPaths_getIndex(loc == 0)')
                error('Invalid group & channel specifications found')
            end
            %NOTE: we copy indexStruct from the GET_INDICES structure
            %thus we've copied the .indices field, still need to define
            %the id field
            indexStruct = GET_INDICES;
            ids = num2cell(loc);
            [indexStruct.id] = deal(ids{:});
        end
        
        %OUTPUT THUS FAR:
        %-------------------------------------
        %indexStruct
        %   .id      -> index of an object to retrieve raw data from
        %   .indices -> start and (stop OR length) of indices (samples) to retrieve           
        
        %Interleaved data check
        %-----------------------------------------------------------------
        isInterleavedSegment = find(logical([segInfo(:).isInterleaved]));
        if ~isempty(isInterleavedSegment)
            if any(ismember(unique([segInfo(isInterleavedSegment).objOrder]),[indexStruct.id]))
                error(['Some channel objects that have subsets of data requested' ...
                    ' are interleaved, this is currently not supported'])
            end
        end
        
        %GETTING READ/SEEK INSTRUCTIONS
        %==================================================================
        %To read subsets of data we are going to form read/seek
        %instructions.  The quickest way of reading subsets of data from a
        %channel or multiple channels is to read what is needed, seek to
        %the next point of needed data, and to read some more.  This next
        %part of code forms these read/seek instructions.  First
        %instructions are generated for each channel, then all of these
        %instructions are interleaved so that we don't need to backtrack.
        %For each read, we specify which channel that particular
        %read belongs to.
        
        %This will occur in four steps
        %1) For each object, reduce subset instructions to clean
        %   start and stop indices -> yLin variable
        %2) Determining which chunk #s contain the start sample and end
        %   sample for each stretch of indices -> bin2d
        %3) For each chunk that is going to be read, specify where to
        %   start, how much to read, and which chunk it is that we are going
        %   to read
        %4) Take these instructions and relate them back to the file
        %Finally, once we have done this for all objects we will sort the
        %reads based on the start read position in the file
        
        numValuesToGetActual = zeros(1,length(metaStruct.numberDataPoints));
        subsetInfo           = cell(1,length(indexStruct));
        for iObject = 1:length(indexStruct)
            
            curEntry = indexStruct(iObject);
            curRaw   = rawDataInfo(curEntry.id);
            
            %STEP 1)
            %==============================================================
            %Setup of x & y
            %--------------------------------------------------------------
            %- x represents the 1st sample of each "chunk"
            %     i.e. which sample of that channel is represented
            %     by that chunk
            %- y are the subsets of the data we are going to grab 
            %- x is the size of the # of chunks that contain that object
            %- y is user input specified
            %- see TDMS_preprocessFile towards end for dataMatrix definition
            %- chunkIndex represents the last valid chunk, since the # of 
            %  rows will be lareger due to preallocation
            nChunksCurObject = curRaw.chunkIndex;
            x = curRaw.dataMatrix(1:nChunksCurObject,2); 
            y = curEntry.indices;
            
            %NOTE: After this point, y(:,2) is now an index
            %regardless of user input intention, this flag calculates index
            %in case the user wants the 2nd input to represent length
            if SUBSET_IS_LENGTH
                y(:,2) = y(:,1) + y(:,2) - 1;
            end
            
            %This little bit of code reduces artificial breaks
            %that the user may have introduced
            %--------------------------------------------------------------
            if size(y,1) > 1
                %Example of y:
                %[5 10;
                %10 15;
                %15 20
                %25 30];
                %
                %Should result in:
                %[5 20;
                %25 30];
                I_fix = find(y(2:end,1) == y(1:end-1,2));
                if ~isempty(I_fix)
                    %Run in reverse to allow multiple merges
                    %Like the example of 5 to 20 above from
                    %5 to 10, 10 to 15, 15 to 20
                    I_fix = I_fix(:)';
                    for iMerge = I_fix(end:-1:1)
                        y(iMerge,2) = y(iMerge+1,2);
                        y(iMerge+1,:) = [];
                    end
                end
            end
            
            nRows = size(y,1);
            yLin  = y';
            yLin  = yLin(:);
            %Now: yLin(1) = start index 1, 
            %     yLin(2) = end index 1
            %     yLin(3) = start index 2, 
            %     yLin(4) = end index 2, etc
            
            
            %ERROR CHECKING
            %-------------------------------------------------------------
            if ~issorted(yLin)
                error('index inputs for GET_INDICES should be sorted')
                %NOTE: for SUBSET_GET this shouldn't be an issue
                %due to the way that the request is passed in
            end
            
            if yLin(1) < 1 || yLin(end) > metaStruct.numberDataPoints(curEntry.id)
                fprintf('ERROR FOR: %s\n',objectPaths{curEntry.id})
                error('indices requested are out of range for current object')
            end
            
            %STEP 2
            %==============================================================
            %Calculation of bin2d, nGrabs, and nChunksPerRow
            %--------------------------------------------------------------
            %NOTE: We add 0.5 so that histc will include an element, as
            %histc normally does a check for  EDGES(k) <= X(i) < EDGES(k+1)
            %NOT,  EDGES(k) <= X(i) <= EDGES(k+1)
            yLin(2:2:end) = yLin(2:2:end) + 0.5;
            [~,bin]       = histc(yLin,x);
            %bin, represents which chunks of data each .indices value
            %belongs to, both for the start (column 1, and the finish
            %column 2), If these two values are not the same, then the
            %subset of data requested in .indices for a particular row
            %spans multiple data chunks in the file
            bin2d = reshape(bin(:),[2 nRows])';
            %bin2d 
            %   - column 1, index of which chunk has the first data
            %     point for each row in indices
            %   - column 2, index of which chunk has the last data
            %     point for each row in indices
            
            nChunksPerRow = (bin2d(:,2) - bin2d(:,1)) + 1; 
            nGrabs        = sum(nChunksPerRow); 
            %NOTE: nGrabs represents the # of reads to perform

            %STEP 3
            %==============================================================
            %Parsing read instructions -> start, length, which chunk
            %--------------------------------------------------------------
            giChunkNumbers   = zeros(nGrabs,1);  %gi -> grab info
            giNumSamplesRead = zeros(nGrabs,1);
            giSampleNumber   = zeros(nGrabs,1);
                        
            curIndex = 0;
            for iSet = 1:nRows 
                if nChunksPerRow(iSet) == 1 
                    %- first and last sample from indices are in the same chunk
                    %- only grab data from 1 chunk
                    curIndex = curIndex + 1;
                    giChunkNumbers(curIndex)   = bin2d(iSet,1);
                    giNumSamplesRead(curIndex) = y(iSet,2) - y(iSet,1) + 1;
                    giSampleNumber(curIndex)   = y(iSet,1);
                else
                    curBin    = bin2d(iSet,1);
                    startSamp = y(iSet,1);
                    lastSamp  = x(curBin+1) - 1;
                    
                    %Notes: 
                    %- we only go up to the last chunk to grab (hence the -1)
                    %- all grab lengths are referenced to the end of the chunk 
                    for iD = 1:nChunksPerRow(iSet)-1
                        curIndex = curIndex + 1;
                        giChunkNumbers(curIndex)   = curBin;
                        giNumSamplesRead(curIndex) = lastSamp - startSamp + 1;
                        giSampleNumber(curIndex)   = startSamp;
                        
                        startSamp = lastSamp + 1;
                        curBin    = curBin + 1;
                        lastSamp  = x(curBin+1) - 1;
                    end
                    %- the last one is referenced to the last sample to grab
                    curIndex = curIndex + 1;
                    giChunkNumbers(curIndex,1) = curBin;
                    giNumSamplesRead(curIndex) = y(iSet,2) - startSamp + 1;
                    giSampleNumber(curIndex)   = startSamp;
                end
            end
            
            %STEP 4
            %==============================================================
            %COLUMN DEFINITIONS
            %==================
            %1) positon to start read, unless a string, then this points
            %   to the start of the chunk, where we get instructions
            %   on how long each string is
            %2) # of samples to grab in that chunk (i.e. read length)
            %3) index of first read sample, relative to chunk start
            %   this is needed for strings :/
            %4) channelID, index of channel for later use in data saving
            %5) seg #, which segment the chunk belongs to, only used
            %   for strings :/
            
            grabInfo = zeros(nGrabs,5);
            nSampsAboveStart = giSampleNumber - x(giChunkNumbers);
            fileStartPos     = curRaw.dataMatrix(1:nChunksCurObject,1);
            if curRaw.dataType == 32
                grabInfo(:,1) = fileStartPos(giChunkNumbers);
            else
                grabInfo(:,1) = fileStartPos(giChunkNumbers) + nSampsAboveStart.*TDMS_getDataSize(curRaw.dataType);
            end
            
            grabInfo(:,2) = giNumSamplesRead;
            grabInfo(:,3) = nSampsAboveStart + 1;
            grabInfo(:,4) = curEntry.id; %This will be needed with multiple objects
            
            segNumbers       = curRaw.dataMatrix(1:nChunksCurObject,3);
            grabInfo(:,5)    = segNumbers(giChunkNumbers); %needed for strings :/
            
            %Put temporarily into a cell array
            subsetInfo{iObject} = grabInfo;
            
            %This gets used in data initialization & error checking
            numValuesToGetActual(curEntry.id) = sum(giNumSamplesRead);
        end %END OF LOOP OVER EACH OBJECT
        %==================================================================
        
        
        if length(indexStruct) == 1
            subsetInfo = grabInfo;
        else
            subsetInfo = cat(1,subsetInfo{:});
            %Resort rows by read order, so that we minimize the distance of fseeks
            [~,I] = sort(subsetInfo(:,1));
            subsetInfo = subsetInfo(I,:);
        end
    end
    
    optionStruct = struct('keepDataArray',keepDataArray,'useSubset',useSubset,...
        'subsetInfo',subsetInfo,'numValuesToGetActual',numValuesToGetActual);
    varargout{1} = optionStruct;
end