function output = TDMS_dataToGroupChanStruct_v4(inputStruct,varargin)
%TDMS_dataToGroupChanStruct_v4  
%
%   NOTE: This reproduces the data structure for the old lab functionality 
%   of the function getStructTDM
%   
%   See Also: TDMS_genvarname2, TDMS_readTDMSFile

REPLACE_STR = '_';
PREPEND_STR = 'p_';
PREPEND_GRP_STRING = 'g_';
PREPEND_CHAN_STRING = 'c_';
ALWAYS_PREPEND = false;
PROP_NAME = 'Props';

propNames    = inputStruct.propNames;
propValues   = inputStruct.propValues;
groupIndices = inputStruct.groupIndices;
groupNames   = inputStruct.groupNames;
chanIndices  = inputStruct.chanIndices;
chanNames    = inputStruct.chanNames;
rootIndex    = inputStruct.rootIndex;
data         = inputStruct.data;

groupNames = cellfun(@(x) TDMS_genvarname2(x,REPLACE_STR,PREPEND_GRP_STRING,ALWAYS_PREPEND),groupNames,'UniformOutput',false);



rootPropStruct = propsToStruct(propNames{rootIndex},propValues{rootIndex},REPLACE_STR,PREPEND_STR,ALWAYS_PREPEND);
output = struct(PROP_NAME,rootPropStruct);

for iGroup = 1:length(groupIndices)
    curGroupIndex = groupIndices(iGroup);
    curChanIndices = chanIndices{iGroup};
    curChanNames   = chanNames{iGroup};
    
    curChanNames = cellfun(@(x) TDMS_genvarname2(x,REPLACE_STR,PREPEND_CHAN_STRING,ALWAYS_PREPEND),curChanNames,'UniformOutput',false);
    
    groupPropStruct = propsToStruct(propNames{curGroupIndex},propValues{curGroupIndex},REPLACE_STR,PREPEND_STR,ALWAYS_PREPEND);
    groupStruct = struct('name',groupNames(iGroup),PROP_NAME,groupPropStruct);
    for iChan = 1:length(curChanIndices)
        curChanIndex = curChanIndices(iChan);
        chanPropStruct = propsToStruct(propNames{curChanIndex},propValues{curChanIndex},REPLACE_STR,PREPEND_STR,ALWAYS_PREPEND);
        chanStruct =  struct('name',curChanNames{iChan},PROP_NAME,chanPropStruct,...
            'data',[]);
        chanStruct.data = data{curChanIndex};
        %NOTE: I had a case statement in case the data type was a string,
        %which would change the interpretation of the struct input
        if strcmp(chanStruct.name,'name')
            error('The variable "name" for a channel is off limits, need a different conversion wrapper (probably 2)')
        end
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

