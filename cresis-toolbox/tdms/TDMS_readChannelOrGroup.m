function [channelData,channelNames] = TDMS_readChannelOrGroup(TDMS_filename,group,channel,subset)
%TDMS_readChannelOrGroup  Return data from a specific channel or group of channels
%
%   Original Author: Preston K. Manwaring
%   Last modified  : 4/29/2011 JAH
%
%   FORMAT 0:   First call to examine file contents
%   =======================================================================
%   [tempStruct] = TDMS_readChannelOrGroup(TDMS_filename,1)
%   
%   
%   NOTE: For FORMAT 1 & FORMAT 2, TDMS_filename can be replaced with
%   the output from the FORMAT 0 call, tempStruct, to save on processing 
%   time
%   
%   FORMAT 1:   To retrieve all channels in a group
%   =======================================================================
%   [channelData,channelNames] = 
%               TDMS_readChannelOrGroup(TDMS_filename,groupName,[],*subset)
%   
%       ALTERNATIVELY ...
%
%   [...] = TDMS_readChannelOrGroup(tempStruct,groupName,[],*subset)
%
%
%   FORMAT 2:   To retrieve a specific channel from a group
%   =======================================================================
%   channelData = TDMS_readChannelOrGroup(TDMS_filename,groupName,channelName,*subset)
%
%
%   Description: This file will read a channel out of a TDMS file and write
%   the data to a local variable. It is a wrapper for TDMS_readTDMSfile.
%
%   INPUTS:
%   =======================================================================
%   tempStruct    : output from previous call to this function, see
%                   definition in OUTPUTS
%   TDMS_filename : The name of the TDMS filename as a string 
%                   (e.g.'C:\foo.tdms')
%   groupName     : The group name in which the channel data reside 
%                   (e.g. 'data Group')

%
%   OPTIONAL INPUTS
%   =======================================================================
%   channelName   : (default []), if empty, retrieves all channels in a
%                   group, otherwise specify the channel name to retrieve
%                   (e.g. 'channel A')
%   subset        : A subset of the raw data to retrieve, format is:
%                   [starting sample  , # of samples to retrieve]
%                   e.g. [100 50] would start retrieving data at the 100th
%                   index and would retrieve 50 samples, 100 through 149
%
%   OUTPUTS:
%   =======================================================================
%   tempStruct    : (structure array with fields)
%       .tdmsStruct    - main output of TDMS_readTDMSfile with no raw data 
%                        retrieved
%       .metaStruct    - second output of TDMS_readTDMSfile
%       .TDMS_filename - the filename that was passed in, this is used if 
%                        this structure gets passed in as an input
%   channelData   : If requesting a single channel, channelData is the
%                   data from that channel.  If requesting a group,
%                   channelData is a cell array, where each index holds the
%                   data for a specific channel
%   channelNames  : (cell array of strings) Name of each channel, used to
%                   match channelData to a channel name
%              
%   EXAMPLES
%   ================================================================
%   1) Read 'myChannel' from 'myGroup', return all raw data
%   channelData = TDMS_readChannelOrGroup('C:\myFile.tdms','myGroup','myChannel')
%
%   2) Retrieve all channels from myGroup, returning only the first 100 values
%   channelData = TDMS_readChannelOrGroup('C:\myFile.tdms','myGroup',[],[1 100])
%
%   See also:
%   TDMS_readTDMSFile

if exist('subset','var')
    SUBSET_GET = subset;
else
    SUBSET_GET = [];
end

if ~exist('channel','var')
    channel = [];
end


%This is a bit messy but it works ...
if isstruct(TDMS_filename)
    %Unpack from previous call
    tempStruct = TDMS_filename;
    
    tdmsStruct    = tempStruct.tdmsStruct;
    metaStruct    = tempStruct.metaStruct;
    TDMS_filename = tempStruct.TDMS_filename;
else
    % get the TDMS file name properties before trying to load a big file
    [tdmsStruct,metaStruct] = TDMS_readTDMSFile(TDMS_filename,'GET_DATA_OPTION','getNone');
    if isnumeric(group)
        channelData = struct('tdmsStruct',tdmsStruct,'metaStruct',metaStruct,'TDMS_filename',TDMS_filename);
        return
    end
end


% check to see if the group name given is valid. If not, give alternatives.
I = find(strcmpi(group,tdmsStruct.groupNames),1);
if isempty(I)
    fprintf(2,'Available group names are:\n');
    disp(char(tdmsStruct.groupNames));
    error('Specified group name: %s, not found in TDMS file.',group);
end

% check to see if the channel name is valid. If not, give alternatives.
if ~isempty(channel)
    J = find(strcmpi(channel,tdmsStruct.chanNames{I}),1);
    if isempty(J)
        fprintf(2,'Available chan names are:\n');
        disp(char(tdmsStruct.chanNames{I}));
        disp('')
        error('Specified channel name: %s, not found in TDMS file.',channel);
    end
end

% Populate request structure
if ~isempty(channel)
    indexOriginalPath      = tdmsStruct.chanIndices{I}(J);
    getStruct = struct(...
        'fullPathsKeep',tdmsStruct.objectPathsOrig{indexOriginalPath});
else
    getStruct = struct('groupsKeep',group);
end

%Retrieve data
tdmsStruct = TDMS_readTDMSFile(TDMS_filename,...
    'META_STRUCT',metaStruct, ...
    'GET_DATA_OPTION','getSubset',...
    'OBJECTS_GET',getStruct, ...
    'SUBSET_GET',SUBSET_GET);

%OUTPUT SETUP
%==============================
if ~isempty(channel)
    channelData = tdmsStruct.data{indexOriginalPath};
    channelNames = {channel};
else
    nChans = length(tdmsStruct.chanNames{I});
    channelData = cell(1,nChans);
    chanIndices = tdmsStruct.chanIndices{I};
    for iChan = 1:nChans
        channelData(iChan) = tdmsStruct.data(chanIndices(iChan));
    end
    channelNames = tdmsStruct.chanNames{I};
end

end

