function [outputStruct,navigator] = TDMS_dataToGroupChanStruct_v2(inputStruct,varargin)
%TDMS_dataToGroupChanStruct_v2  Converts output into tiered structure
%
%   [outputStruct,navigator] = TDMS_dataToGroupChanStruct_v2(inputStruct)
%
%   This function tiers the data into groups & channels but leaves the
%   channel & group names as properties of the array, as opposed to making
%   each group and channel its own field name.
%
%   For example for groups => {'servers' 'computers'}
%   
%   Instead of indexing a.servers & a.computers you have to do
%   s.group(1).id => 'servers'
%   s.group(2).id => 'computers'
%
%   This may be desirable for cases in which the names are not appropriate
%   variable names.  The navigator output is provided for going from a name
%   to an index into the structure.
%
%   INPUTS
%   ============================
%   inputStruct : see TDMS_readTDMSfile
%
%   OUTPUTS
%   ============================
%   outputStruct :
%       .propNames  - cell array of strings (for file)
%       .propValues - cell array, for each string, holds value
%       .groups (structure array)
%           .id         - channel name
%           .propNames  - ""
%           .propValues - ""
%           .chans      - structure array
%               .id           - channel name
%               .propNames    - ""
%               .propValues   - ""
%               .dataTypeName - class name
%               .dataType     - Labview enum to indicate type
%               .data         - data array
%
%   navigator :
%       .groups
%           .names - names of all the groups
%           .indices - indices in which they reside in outputStruct.groups
%       .chans
%           .chanNames     - cell array of strinsg, name of the channel
%           .groupNames    - cell array of strings, group name of the
%                            channel
%           .fullObjNames  - cell array of strings, the full name of the
%                            channel object (as tracked in Labview)
%           .groupIndices
%           .chanIndices   - together with groupIndices allows one to grab
%                            a particular channel object
%
%   See also: TDMS_readTDMSFile

%EXAMPLE OF OUTPUT FORMAT
%==========================================================================
% inputStruct =
%           rootIndex: 1
%        groupIndices: 2
%          groupNames: {'Data'}
%         chanIndices: {[3 4 5]}
%           chanNames: {{1x3 cell}}
%                data: {[]  []  [1x9785383 single]  [1x9785383 single]  [1x9785383 single]}
%           propNames: {{1x0 cell}  {1x0 cell}  {1x1 cell}  {1x1 cell}  {1x1 cell}}
%          propValues: {{1x0 cell}  {1x0 cell}  {1x1 cell}  {1x1 cell}  {1x1 cell}}
%     objectPathsOrig: {'/'  '/'Data''  '/'Data'/'Trigger''  '/'Data'/'Chan1''  [1x15 char]}
%            dataType: [0 0 9 9 9]
%        dataTypeName: {'void'  'void'  'single'  'single'  'single'}


%Root index building
%=========================================
rootIndex = inputStruct.rootIndex;
outputStruct = struct('propNames',[],'propValues',[],'groups',[]);
outputStruct.propNames  = inputStruct.propNames{rootIndex};
outputStruct.propValues = inputStruct.propValues{rootIndex};

%Groups & channels
%=========================================
nChans = sum(cellfun('length',inputStruct.chanIndices));

%These are for the navigator
chanIndices = zeros(1,nChans);
groupIndices = zeros(1,nChans);
chanNames   = cell(1,nChans);
groupNames  = cell(1,nChans);
fullObjNames = cell(1,nChans);
chanCount   = 0;

groups = struct('id',inputStruct.groupNames,'chans',[],'propNames',[],'propValues',[]);
for iGroup = 1:length(groups)
    groupIndex = inputStruct.groupIndices(iGroup);
    groups(iGroup).propNames  = inputStruct.propNames{groupIndex};
    groups(iGroup).propValues = inputStruct.propValues{groupIndex};
    
    %Channel entries
    %-----------------------------------
    chans =  struct('id',inputStruct.chanNames{iGroup},'dataTypeName',[],'dataType',[],'propNames',[],'propValues',[],'data',[]);
    for iChan = 1:length(chans)
        chanIndex = inputStruct.chanIndices{iGroup}(iChan);
        chanCount = chanCount + 1;
        chanIndices(chanCount) = iChan;
        groupIndices(chanCount) = iGroup;
        
        chanNames{chanCount}   = chans(iChan).id;
        groupNames{chanCount}  = groups(iGroup).id;
        fullObjNames{chanCount} = ['/''' groupNames{chanCount} '''/''' chanNames{chanCount} ''''];
        
        chans(iChan).propNames    = inputStruct.propNames{chanIndex};
        chans(iChan).propValues   = inputStruct.propValues{chanIndex};
        chans(iChan).data         = inputStruct.data{chanIndex};
        chans(iChan).dataType     = inputStruct.dataType(chanIndex);
        chans(iChan).dataTypeName = inputStruct.dataTypeName{chanIndex};
    end
    groups(iGroup).chans = chans;
end
outputStruct.groups = groups;

%Build the navigator structure
%====================================================================
%   navigator
%   ===========================
%   - groups
%       - name array
%       - group index
%   - chan name
%       - chan name
%       - group name
%       - group index
%       - chan index


groupsNav = struct('names',inputStruct.groupNames,'indices',num2cell(1:length(groups)));
chansNav  = struct('chanNames',{chanNames},'groupNames',{groupNames},'fullObjNames',...
    {fullObjNames},'groupIndices',groupIndices,'chanIndices',chanIndices);
navigator = struct('groups',groupsNav,'chans',chansNav);

