function [status,message] = opsReleaseLayerGroup(sys,param)
%
% [status,message] = opsReleaseLayerGroup(sys,param)
%
% Sets the public status of a layer group to true.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.lyr_group_name = string
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   message: status message
%
% Author: Kyle W. Purdon

% CONSTRUCT THE JSON STRUCTURE
jsonSruct = struct('properties',param.properties);

% CONVERT THE JSON STRUCUTRE TO A JSON STRING
try
  jsonStr = tojson(jsonSruct);
catch ME
  jsonStr = savejson('',jsonSruct,'FloatFormat','%2.10f','NaN','null');
end

% SEND THE COMMAND TO THE SERVER
opsCmd;
if gOps.profileCmd
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'profile'),gOps.dbUser,gOps.dbPswd,...
    'Post',{'app' sys 'data' jsonStr 'view' 'releaseLayerGroup'});
else
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'release/layer/group'),gOps.dbUser,gOps.dbPswd,...
    'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,message] = jsonResponseDecode(jsonResponse);

end