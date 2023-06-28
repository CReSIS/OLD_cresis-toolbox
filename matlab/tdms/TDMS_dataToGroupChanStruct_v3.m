function output = TDMS_dataToGroupChanStruct_v3(inputStruct,varargin)
%TDMS_dataToGroupChanStruct_v3  
%
%   translates objects AND properties to valid field names
%
%   Similar to TDMS_dataToGroupChanStruct_v1, but also translates
%   properties
%
%   output = TDMS_dataToGroupChanStruct_v3(inputStruct)
%
%   OPTIONAL INPUTS (property/value pairs)
%   ============================================
%   REPLACE_STR = default '_', what to replace invalid characters with
%   PREPEND_STR = default 'v' (for variable), gets prepended if first
%               character is invalid (must be alpha, not numeric or _)
%   ALWAYS_PREPEND = default false, if true always prepends the prepend string,
%               regardless of whether or not it is needed
%
%   
%
%   See Also: TDMS_genvarname2

REPLACE_STR = '_';
PREPEND_STR = 'v';
ALWAYS_PREPEND = false;

propNames    = inputStruct.propNames;
propValues   = inputStruct.propValues;
groupIndices = inputStruct.groupIndices;
groupNames   = inputStruct.groupNames;
chanIndices  = inputStruct.chanIndices;
chanNames    = inputStruct.chanNames;
rootIndex    = inputStruct.rootIndex;
data         = inputStruct.data;

rootPropStruct = propsToStruct(propNames{rootIndex},propValues{rootIndex},REPLACE_STR,PREPEND_STR,ALWAYS_PREPEND);
output = struct('props',rootPropStruct);

for iGroup = 1:length(groupIndices)
    curGroupIndex = groupIndices(iGroup);
    curChanIndices = chanIndices{iGroup};
    curChanNames   = chanNames{iGroup};
    
    groupPropStruct = propsToStruct(propNames{curGroupIndex},propValues{curGroupIndex},REPLACE_STR,PREPEND_STR,ALWAYS_PREPEND);
    groupStruct = struct('name',groupNames(iGroup),'props',groupPropStruct);
    for iChan = 1:length(curChanIndices)
        curChanIndex = curChanIndices(iChan);
        chanPropStruct = propsToStruct(propNames{curChanIndex},propValues{curChanIndex},REPLACE_STR,PREPEND_STR,ALWAYS_PREPEND);
        chanStruct =  struct('name',curChanNames{iChan},'props',chanPropStruct,...
            'data',[]);
        chanStruct.data = data{curChanIndex};
        %NOTE: I had a case statement in case the data type was a string,
        %which would change the interpretation of the struct input
        groupStruct.(TDMS_genvarname2(chanStruct.name,...
            REPLACE_STR,PREPEND_STR,ALWAYS_PREPEND)) = chanStruct;
    end
    output.(TDMS_genvarname2(groupStruct.name,...
        REPLACE_STR,PREPEND_STR,ALWAYS_PREPEND)) = groupStruct;
end
end

function propStruct = propsToStruct(names,values,REPLACE_STR,PREPEND_STR,ALWAYS_PREPEND)
propStruct = struct([]);
for iProp = 1:length(names)
   propStruct(1).(TDMS_genvarname2(names{iProp},REPLACE_STR,PREPEND_STR,ALWAYS_PREPEND)) = values{iProp}; 
end
end
