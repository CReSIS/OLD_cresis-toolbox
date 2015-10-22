function [status,data] = opsGetLayers(sys)
%
% [status,data] = opsGetLayers(sys)
%
% Retrieves all layers for a given system with normal status from the database.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: structure with fields (or error message)
%       properties.lyr_id = double array
%       properties.lyr_name = cell array
%       properties.lyr_group_name = cell array
%
% Author: Kyle W. Purdon, Trey Stafford

global gRadar;

% CONSTRUCT THE JSON STRUCTURE
opsCmd;
param.properties.mat = true;
opsAuth = load(fullfile(gRadar.tmp_path,'ops.mat'));
param.properties.userName = opsAuth.userName;
param.properties.isAuthenticated = opsAuth.isAuthenticated;
jsonStruct = struct('properties',param.properties);

% CONVERT THE JSON STRUCTURE TO A JSON STRING
try
  jsonStr = tojson(jsonStruct);
catch ME
  jsonStr = savejson('',jsonStruct,'FloatFormat','%2.10f','NaN','null');
end

% SEND THE COMMAND TO THE SERVER
if gOps.profileCmd
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'profile'),'','',...
    'Post',{'app' sys 'data' jsonStr 'view' 'getLayers'});
else
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'get/layers'),'','',...
    'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,decodedJson] = jsonResponseDecode(jsonResponse);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
data.properties.lyr_name = decodedJson.lyr_name;
data.properties.lyr_id = double(cell2mat(decodedJson.lyr_id));
data.properties.lyr_group_name = decodedJson.lyr_group_name;

end