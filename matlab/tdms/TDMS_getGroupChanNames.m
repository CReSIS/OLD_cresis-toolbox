function metaStruct = TDMS_getGroupChanNames(metaStruct)
%TDMS_getGroupChanNames  Small function to get group/channel names from paths
%
%   metaStruct = TDMS_getGroupChanNames(metaStruct)
%
%   INPUT
%   =======================================================================
%   metaStruct (structure) containing field:
%   	.objectNameList : (cellstr), name of each object (unmodified from
%                          file) of format => /'<group>'/'<chan>'
%
%   OUTPUT
%   =======================================================================
%   metaStruct (structure) with added fields:
%       .groupNames : (cellstr), same length as objectNameList, group name
%                      for each input object, root will be empty
%       .chanNames  : (cellstr), "     ", channel name for each object,
%                      root and group objects will be empty
%       .isChan     : (logical)
%
%   EXAMPLE OUTPUT (of relevant fields)
%   =======================================================================
%     objectNameList: {1x68 cell}
%         groupNames: {1x68 cell}
%          chanNames: {1x68 cell}
%             isChan: [1x68 logical]
%
%   See Also: TDMS_readTDMSFile

%PARSING RULES:
%=============================================================
%A path follows the following format:
%channel objects /'<group>'/'<chan>'
%groups objects  /'<group>'
%root            '/'

%ACTUAL SPECIFICATION =====================================================
%http://zone.ni.com/devzone/cda/tut/p/id/5696
%
%Every TDMS object is uniquely identified by a path. Each path is a string
%including the name of the object and the name of its owner in the TDMS
%hierarchy, separated by /. Each name is enclosed by the ' ' symbols. Any'
%symbol within an object name is replaced with two ' symbols. The following
%table illustrates path formatting examples for each type of TDMS object:
%
%It might be possible to get by the current filter ..., CONSIDER REWRITING
%/'myGroup''/''name'/' <- 

objectPaths = metaStruct.objectNameList;

pat1 = '/''(?<groupNames>.*?)''/''(?<chanNames>.*?)''$'; %Channel Object
pat2 = '/''(?<groupNames>.*?)''$(?<chanNames>)';  %Group Object
pat3 = '/(?<groupNames>)(?<chanNames>)';  %Root object
pat = [pat1 '|' pat2 '|' pat3];
temp = regexp(objectPaths,pat,'names');

metaStruct.groupNames = cellfun(@(x) x.groupNames,temp,'UniformOutput',false);
metaStruct.chanNames  = cellfun(@(x) x.chanNames,temp,'UniformOutput',false);
metaStruct.isChan     = ~cellfun('isempty',metaStruct.chanNames);



