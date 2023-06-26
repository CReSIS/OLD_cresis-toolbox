function [finalOutput,metaStruct] = TDMS_readTDMSFile(tdmsFileName,varargin)
%TDMS_readTDMSFile  Reads TDMS file and does minimal processing to obtain output
%
%   [finalOutput,metaStruct] = TDMS_readTDMSFile(tdmsFileName)
%
%   [...] = TDMS_readTDMSFile(tdmsFileName,'Property1',PropertyValue1,...)
%   allows for specification of additional properties.
%
%   [...] = TDMS_readTDMSFile(tdmsIndexFileName,...)  Allows you to pass 
%   in the tdms_index file for debugging
%
%   This is the main file for reading TDMS files.  It reads the file and
%   does minimal processing on the output.  
%
%   For wrappers of this function (RECOMMENDED):
%   See Also:
%       TDMS_getStruct
%       TDMS_dataToGroupChanStruct_v1,
%       TDMS_dataToGroupChanStruct_v2,
%       TDMS_dataToGroupChanStruct_v3,
%       TDMS_dataToGroupChanStruct_v4
%       TDMS_readChannelOrGroup
%
%   RETRIEVING PARTIAL DATA
%   =======================================================================
%   This code has the ability to retrieve a subset of the data from the
%   tdms file.  For documentation on this see: TDMS_retrievingSubsets
%
%   OUTPUTS: 
%   NOTE: Wrappers exist to change this to a preferred format.  See above.
%   =======================================================================
%   finalOutput : (struct)
%   ---------------------------------
%       .rootIndex    - index of the root
%     	.groupIndices - (numeric array) indices of all group objects
%      	.chanIndices  - (cell array of numeric arrays)
%                       ex. {[1 2]   [3 4 5]   [6 7 8]}
%                       Each index in the cell array corresponds to the
%                       same index in groupNames and groupIndices
%                       At that index is an array of indices for the
%                       channel objects that belong to that group
%      	.chanNames    - (cell array of cell array of strings)
%                       ex. {{'1' '2'} {'3' '4' '5'}}
%                       similar setup to chanIndices
%     	.groupNames   - (cell array of strings) names of each group object
%
%           For the following properties, the indices referred to above , 
%           groupIndices and chanIndices, index into these arrays:
%      	.data         - (cell array, # elements = # objects)
%     	.propNames    - (cell array of cell arrays)
%       .propValues   - (cell arrays of cell arrays)
%       .dataType     - (num. array) Labview dataType enumeration
%       .dataTypeName - (cell array of strings)
%       .objectPathsOrig - (cell array of strings) This is thrown in for 
%                          reference and is the full name of the object as
%                          tracked in the TDMS file
%       .numberDataPointsRaw - The # of data points for each object that
%                              are available in the TDMS file. This is not
%                              filtered, thus it may not match the # of
%                              data points returned in .data
%
%   metaStruct : (struct)
%   ---------------------------------
%       This output is essentially a dumping ground for a lot of the
%   temporary variables in the code and may change between versions.
%       .version   - indicates which version of the code was used
%                    to create this structure
%       .fileName  - relative name (no path) of the tdms file on which
%                    this structure is based
%       .numberDataPoints - # of data points availabe for each channel
%       .chanNames - channel name for each object, empty for channels &
%                    root
%       .groupNames - group name for each object
%   OTHER VARIABLES IN metaStruct MAY CHANGE ...
%   
%
%   OPTIONAL PARAMETERS: pass in as property/value pairs, case insensitive
%   =======================================================================
%   UTC_DIFF        : (default -5) Conversion of UTC Timestamp to local
%                     time, if desired this can be set to 0 and no 
%                     conversion will occur
%                           
%                     Eastern Time: -5, Central Time: -6
%
%   USE_INDEX       : (default true) If true, reads the tdms_index file
%                     if it is available.  This may offer some speedup in
%                     reading the file.
%
%   MAX_NUM_OBJECTS : (default 100) This is currently used in preallocation
%                      of the # of objects.  On overflow the same # is used
%                      to grow by, i.e. on object 101, grow to 200
%
%   MAX_NUM_PROPS   = (default 20), used for initializing property cell
%                      arrays for each object
%
%   N_SEGS_GUESS    = (default 25000), estimated # of segments for
%                      preallocation purposes
%
%   N_SEGS_INC      = (default 25000), how much to expand upon segment
%                      overflow, i.e. if more than N_SEGS_GUESS segments
%
%   DATE_STR_FORMAT = (default 'dd-mmm-yyyy HH:MM:SS:FFF'), how to process
%                     the timestamp properties, timestamp arrays (i.e.
%                     data) are left in the same format as Matlab's now
%                     command, i.e. a numeric value that you can call the
%                     function datestr() on to return a value
%
%   INIT_CHUNK_SIZE = (default 1000), the expected # of chunks of data,
%                     this probably shouldn't be changed ...
%
%   Additional properties are defined in TDMS_retrievingSubsets



%FUNCTIONS CALLED IN THIS FILE:
%==============================================
%       TDMS_handleGetDataOption
%       TDMS_preprocessFile
%       TDMS_getGroupChanNames
%       TDMS_processLeadIn
%       TDMS_readFileHelper_v2
%       TDMS_readFileHelper_v1


if nargin == 0
    error('Input to %s requires a filename',mfilename)
end

if ~exist(tdmsFileName,'file')
    error('The file specified for tdms reading doesn''t exist, FILENAME: %s',tdmsFileName)
end

%Added to allow for tdms_index reading
[~,~,fileExt] = fileparts(tdmsFileName);

%DON'T CHANGE ME
%===========================================
temp = [who; {'temp'}];
STRING_ENCODING = 'UTF-8';
SECONDS_IN_DAY  = 86400;  %60s * 60 min * 24 hours
CONV_FACTOR     = 695422; %datenum('01-Jan-1904')
ROOT_PATH       = '/';
ROOT_GROUP_NAME = '';
TDMS_INDEX_EXT  = '.tdms_index';
MACHINE_FORMAT  = 'ieee-le'; %NOTE: Eventually this could be passed into
%fread for support of different endianess
CURRENT_VERSION = 2.5;
constantParams = setdiff(who,temp);  %NOTE: Below I place these into the 
%constants variable structure but I don't allow assignment from inputs


%PROPERTIES THAT MIGHT NEED TO BE CHANGED
%==========================================================================
temp = [who; {'temp'}];
%The following get passed down into TDMS_preprocessFile
UTC_DIFF        = -5;   %#ok<*NASGU> %This will report times in Eastern Standard Time
USE_INDEX       = true; %Parse meta data from index if present
MAX_NUM_OBJECTS = 100;  %Best guess as to max # of objects ever encountered
MAX_NUM_PROPS   = 20;   %"                        " properties "          "
N_SEGS_GUESS    = 25000;%same, but for segments
N_SEGS_INC      = 25000;%if off, how large to grow estimate by when resizing
DATE_STR_FORMAT = 'dd-mmm-yyyy HH:MM:SS:FFF'; %default conversion of property value
                  %NOTE: data values are not converted, but run datestr()
                  %on them if you wish to see their value as a string
INIT_CHUNK_SIZE = 1000; %a chunk is a subset of a segment, the # of chunks 
                        %varies on an object by object basis, only applies
                        %to objects with raw data
                                  
DEBUG           = false; %Internal use

%NOTE: INDEX_DEBUG indicates that we'll only attempt to read the index
%file, no raw data will be read
if strcmp(fileExt,'.tdms_index')
    INDEX_DEBUG     = true;
else
    INDEX_DEBUG     = false;
end

OBJECTS_GET      = struct([]);
OBJECTS_IGNORE   = struct([]);
SUBSET_GET       = [];
GET_INDICES      = [];  %Needs fields 'group' 'channel' 'indices'
SUBSET_IS_LENGTH = true;
GET_DATA_OPTION  = 'getAll';
META_STRUCT      = []; %Can be used to reduce processing time
                      %This is an output from this file 
defaultOPTS  = setdiff(who,temp);  %Anything in defaultOPTS can be overwritten
%below if the name matches a property name from the input

%Assigment of default values, as well as paramsStruct population
%--------------------------------------------------------------------------
paramsStruct = struct;
if nargin > 1
    defNames = varargin(1:2:end);
else
    defNames = {};
end

for iDefault = 1:length(defaultOPTS)
    curVariable = defaultOPTS{iDefault};
    I = find(strcmpi(curVariable,defNames),1);
    if ~isempty(I)
        eval([curVariable '= varargin{2*I};'])
    end
    paramsStruct.(curVariable) = eval(curVariable);
end

for iConst = 1:length(constantParams)
    curVariable = constantParams{iConst};
    paramsStruct.(curVariable) = eval(curVariable);
end
%paramsStruct now contains all constants, as well as their values

%Data option checking
%--------------------------------------------------------------------------
TDMS_handleGetDataOption('check',paramsStruct)

%File Opening & meta data processing
%--------------------------------------------------------------------------
if ~INDEX_DEBUG
    fid = fopen(tdmsFileName,'r',MACHINE_FORMAT,STRING_ENCODING);
else
    fid = []; %Don't open the .tdms file, .tdms_index only
end

metaStruct   = TDMS_preprocessFile(fid,tdmsFileName,paramsStruct);
metaStruct   = TDMS_getGroupChanNames(metaStruct); %Processing of channel & group names
optionStruct = TDMS_handleGetDataOption('getArray',paramsStruct,metaStruct); %Option handling


%Propery handling
%--------------------------------------------------------------------------
rawDataInfo = metaStruct.rawDataInfo;
numObjects  = length(rawDataInfo);
propNames   = cell(1,numObjects);
propValues  = cell(1,numObjects);
for iObject = 1:numObjects
    propNames{iObject}  = rawDataInfo(iObject).propNames;
    propValues{iObject} = rawDataInfo(iObject).propValues;
end

%Getting the data
%--------------------------------------------------------------------------
if strcmpi(GET_DATA_OPTION,'getNone') || INDEX_DEBUG
    data = cell(1,length(metaStruct.rawDataInfo));
    if ~INDEX_DEBUG
        fclose(fid);
    end
elseif optionStruct.useSubset
    data = TDMS_readFileHelper_v2(fid,optionStruct,metaStruct,paramsStruct);
else
    data = TDMS_readFileHelper_v1(fid,optionStruct,metaStruct,paramsStruct);
end

%INITIALIZATION OF GENERIC OUTPUT
%==========================================================================
%use groupInfo = TDMSp_getGroupInfo(groupNames,chanNames,isChan) instead

[uGroups,~,Igroup] = unique(metaStruct.groupNames);
%NOTE: Igroup tells us which of the unique group names each index belongs to

rootIndex = find(cellfun(@(x) strcmp(x,ROOT_PATH),metaStruct.objectNameList));

groupIndices     = zeros(1,length(uGroups)-1);
chanIndices      = cell(1,length(groupIndices));
groupNamesOutput = cell(1,length(groupIndices));
chanNamesOutput  = cell(1,length(groupIndices));

curGroupCount = 0;
for iGroup = 1:length(uGroups)
    if ~strcmp(uGroups(iGroup),ROOT_GROUP_NAME)
        curGroupCount = curGroupCount + 1;
        
        %Find all those with matching group name
        inGroup = find(Igroup == iGroup);
        
        %Parse out the actual group object
        isGroupObject = ~metaStruct.isChan(inGroup);
        I_groupObject = inGroup(isGroupObject);
        groupIndices(curGroupCount)     = I_groupObject;
        groupNamesOutput(curGroupCount) = uGroups(iGroup);
        
        %Parse out the channels for that group
        inGroup(isGroupObject)     = [];
        chanIndices{curGroupCount} = inGroup;
        chanNamesOutput{curGroupCount} = metaStruct.chanNames(inGroup);
    end
end

finalOutput = struct(...
    'rootIndex',rootIndex,...
    'groupIndices',groupIndices,...
    'groupNames',{groupNamesOutput},...
    'chanIndices',{chanIndices},...
    'chanNames',{chanNamesOutput},...
    'data',{data},...
    'propNames',{propNames},...
    'propValues',{propValues},...
    'objectPathsOrig',{metaStruct.objectNameList},...
    'numberDataPointsRaw',metaStruct.numberDataPoints,...
    'dataType',[rawDataInfo.dataType],...
    'dataTypeName',{arrayfun(@TDMS_getDataTypeName,[rawDataInfo.dataType],'UniformOutput',false)});
end
