function [status,message] = opsLogoutUser()
%
% [status,data] = opsLogoutUser()
%
% Logs out a user from the OPS system.
%
% Input:
%   none
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: string status message
%
% Author: Kyle W. Purdon

global gRadar;

% CONSTRUCT THE JSON STRUCTURE
opsCmd;
param.properties.mat = true;
opsAuth = load(fullfile(gRadar.tmp_path,'ops.mat'));

if ~opsAuth.isAuthenticated
  error('ERROR: NO USER IS LOGGED IN FROM MATLAB');
end

param.properties.userName = opsAuth.userName;

jsonSruct = struct('properties',param.properties);

% CONVERT THE JSON STRUCUTRE TO A JSON STRING
try
  jsonStr = tojson(jsonSruct);
catch ME
  jsonStr = savejson('',jsonSruct,'FloatFormat','%2.10f','NaN','null');
end

% SEND THE COMMAND TO THE SERVER
if gOps.profileCmd
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'profile'),'','',...
    'Post',{'app' 'rds' 'data' jsonStr 'view' 'logoutUser'});
else
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'logout/user'),'','',...
    'Post',{'app' 'rds' 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,message] = jsonResponseDecode(jsonResponse);

% WRITE OPS AUTH FILE AND CLEAR THE PROFILE
if status == 1
  opsAuth.isAuthenticated = false;
  save(fullfile(gRadar.tmp_path,'ops.mat'),'-struct','opsAuth')
  opsProfile = struct();
  save(fullfile(gRadar.tmp_path,'ops.profile.mat'),'-struct','opsProfile')
end
end