function [status,message] = opsLoginUser()
%
% [status,data] = opsLoginUser()
%
% Logs a user into the OPS system.
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

opsAuthFn = fullfile(gRadar.tmp_path,'ops.mat');
if exist(opsAuthFn,'file')
  opsAuth = load(fullfile(gRadar.tmp_path,'ops.mat'));
else
  opsAuth.isAuthenticated = false;
end

if opsAuth.isAuthenticated
  [password, userName] = passwordEntryDialog('enterUserName', true,...
    'DefaultUserName', opsAuth.userName, 'DefaultPassword', ...
    opsAuth.password, 'ValidatePassword', false, 'CheckPasswordLength', false);
else
  [password, userName] = passwordEntryDialog('enterUserName', true,...
    'DefaultUserName', 'anonymous', 'DefaultPassword', ...
    'anonymous','ValidatePassword', false, 'CheckPasswordLength', false);
end

param.properties.userName = userName;
param.properties.password = password;

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
    'Post',{'app' 'rds' 'data' jsonStr 'view' 'loginUser'});
else
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'login/user'),'','',...
    'Post',{'app' 'rds' 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,message] = jsonResponseDecode(jsonResponse);

if ~exist(gRadar.tmp_path,'dir')
  mkdir(gRadar.tmp_path);
end

% WRITE OPS AUTH FILE
if status == 1
  opsAuth.userName = param.properties.userName;
  opsAuth.password = param.properties.password;
  opsAuth.isAuthenticated = true;
  save(fullfile(gRadar.tmp_path,'ops.mat'),'-struct','opsAuth')
  
  % GET AND WRITE OPS PROFILE FILE
  [~,~] = opsGetUserProfileData();
  
else
  opsAuth.userName = param.properties.userName;
  opsAuth.isAuthenticated = false;
  save(fullfile(gRadar.tmp_path,'ops.mat'),'-struct','opsAuth')
end
end