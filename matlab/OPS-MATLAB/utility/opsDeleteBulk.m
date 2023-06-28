function [status,data] = opsDeleteBulk(sys,param)
%
% [status,data] = opsDeleteBulk(sys,param)
%
% Completely removes data from OPS.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.season: (string) season name for given segments ('2011_Greenland_P3')
%     properties.only_layer_points: (boolean) true deletes only layer points
%     properties.segment = cell of segment name(s)
%         {''}: deletes all segments
%         {'20110331_01'}: deletes the given single segment
%         {'20110331_01','20110331_02'}: deletes both given segments
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: Message indicating deletion status
%
% Author: Trey Stafford, Kyle Purdon

% AUTHENTICATE THE USER
[param,~,opsProfile] = opsAuthenticate(param,true);
if ~opsProfile.bulkDeleteData
    error('USER NOT AUTHORIZED TO BULK DELETE DATA');
end

% CONFIRM INPUTS BEFORE PROCEEDING 
if strcmp(param.properties.segment{1},'')
    segments = 'ALL';
else
    segments = sprintf(' %s;',param.properties.segment{:});
end

if param.properties.only_layer_points
    layerPoints = 'True';
else
    layerPoints = 'False';
end
opsCmd;
confirmParams = {sprintf('SERVER: \t %s',gOps.serverUrl),'',sprintf('SEASON NAME: \t %s',param.properties.season),'',...
  sprintf('SYSTEM NAME: \t %s',sys),'',sprintf(' Only Delete Layer Points?: \t %s',layerPoints),'',...
  sprintf(' SEGMENTS: \t %s',segments),''};

confirmButton = questdlg(confirmParams,'Confirm Settings','YES:DELETE','NO:QUIT','NO:QUIT');

switch confirmButton
  case 'NO:QUIT'
    error('PROCESS STOPPED BY USER.');
end

% CONSTRUCT THE JSON STRUCTURE
jsonStruct = struct('properties',param.properties);

% CONVERT THE JSON STRUCTURE TO A JSON STRING
try
    jsonStr = tojson(jsonStruct);
catch ME
    jsonStr = savejson('',jsonStruct,'FloatFormat','%2.10f','NaN','null');
end

% SEND THE COMMAND TO THE SERVER
if gOps.profileCmd
    [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'profile'),gOps.dbUser,gOps.dbPswd,...
        'Post',{'app' sys 'data' jsonStr 'view' 'bulkDelete'});
else
    [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'delete/bulk'),gOps.dbUser,gOps.dbPswd,...
        'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,decodedJson] = jsonResponseDecode(jsonResponse);

% Return the status of the deletion.
if status==1
    data = decodedJson;
    fprintf('%s\n', decodedJson);
else
    data = decodedJson;
    error(decodedJson);
end