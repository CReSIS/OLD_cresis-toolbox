function [status,data] = opsGetSystemInfo()
%
% [status,data] = opsGetSystemInfo()
%
% Retrieves all systems/seasons/locations from the database
%
% Input:
%   none
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: structure with fields (or error message)
%       properties.system = cell of string/s
%       properties.season = cell of string/s
%       properties.location = cell of string/s
%       properties.public = array of boolean/s
%       opsProfile = structure
%
% Author: Kyle W. Purdon

% CONSTRUCT THE JSON STRUCTURE
[param,~,opsProfile] = opsAuthenticate(struct('properties',[]),true);
jsonStruct = struct('properties',param.properties);

% CONVERT THE JSON STRUCTURE TO A JSON STRING
try
  jsonStr = tojson(jsonStruct);
catch ME
  jsonStr = savejson('',jsonStruct,'FloatFormat','%2.10f','NaN','null');
end

% SEND THE COMMAND TO THE SERVER
opsCmd;
if gOps.profileCmd
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'profile'),'','',...
    'Post',{'app' '' 'data' jsonStr 'view' 'getSystemInfo'});
else
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'get/system/info'),'','',...
    'Post',{'app' 'rds' 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,decodedJson] = jsonResponseDecode(jsonResponse);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
data = struct('properties',struct('systems',[],'seasons',[],'locations',[])); %,'public',[]));
for outIdx = 1:length(decodedJson)
  data.properties.systems{end+1} = decodedJson{outIdx}.system;
  data.properties.seasons{end+1} = decodedJson{outIdx}.season;
  data.properties.locations{end+1} = decodedJson{outIdx}.location;
%   data.properties.public(end+1) = decodedJson{outIdx}.public;

end
data.opsProfile = opsProfile;
end