function [status,message] = opsCreateUser()
%
% [status,data] = opsCreateUser(param)
%
% Creates a new user with default permissions in the OPS system.
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
param.properties.mat = true;

[password, userName] = passwordEntryDialog('enterUserName', true, ...
  'DefaultUserName', 'anonymous', 'DefaultPassword', ...
    'anonymous', 'ValidatePassword', false, 'CheckPasswordLength', false, ...
    'WindowName', 'New User');

prompt = {'Email:'};
dlg_title = 'New User';
num_lines = [1, length(dlg_title)+40];
defaults={'anonymous@cresis.ku.edu'};
options.Resize='on';
results = inputdlg(prompt,dlg_title,num_lines,defaults,options);

param.properties.userName = userName;
param.properties.password = password;
param.properties.email = results{1};

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
    'Post',{'app' 'rds' 'data' jsonStr 'view' 'createUser'});
else
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'create/user'),gOps.dbUser,gOps.dbPswd,...
    'Post',{'app' 'rds' 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,message] = jsonResponseDecode(jsonResponse);

% WRITE OPS AUTH FILE
if status == 1
  opsAuth.userName = param.properties.userName;
  opsAuth.password = param.properties.password;
  opsAuth.isAuthenticated = false;
  save(fullfile(gRadar.tmp_path,'ops.mat'),'-struct','opsAuth')
else
  opsAuth.userName = param.properties.userName;
  opsAuth.isAuthenticated = false;
  save(fullfile(gRadar.tmp_path,'ops.mat'),'-struct','opsAuth')
end

end