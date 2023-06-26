function [status,message] = opsGetUserProfileData()
%
% [status,message] = opsGetUserProfileData()
%
% Gets a user profile from the OPS system for a logged in user.
%
% Input:
%   none
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   message: status message
%   data: this is written to the gRadar.tmp_path/opsAuth.mat file
%       properties.rds_season_groups = (string or cell of strings)
%       properties.rds_season_group_ids = (integer or array of integers)
%       properties.rds_layer_groups = (string or cell of strings)
%       properties.rds_layer_group_ids = (integer or array of integers)
%       properties.accum_season_groups = (string or cell of strings)
%       properties.accum_season_group_ids = (integer or array of integers)
%       properties.accum_layer_groups = (string or cell of strings)
%       properties.accum_layer_group_ids = (integer or array of integers)
%       properties.snow_season_groups = (string or cell of strings)
%       properties.snow_season_group_ids = (integer or array of integers)
%       properties.snow_layer_groups = (string or cell of strings)
%       properties.snow_layer_groups_ids = (integer or array of integers)
%       properties.kuband_season_groups = (string or cell of strings)
%       properties.kuband_season_group_ids = (integer or array of integers)
%       properties.kuband_layer_groups = (string or cell of strings)
%       properties.kuband_layer_group_ids = (integer or array of integers)
%       properties.layerGroupRelease = (boolean) can the user release layergroups?
%       properties.seasonRelease = (boolean) can the user release seasons?
%       properties.createData = (boolean) can the user create new data?
%       properties.bulkDeleteData = (boolean) can the user bulk delete data?
%       properties.isRoot = (boolean) is the user a root user?
%
% Author: Kyle W. Purdon

global gRadar;

% CONSTRUCT THE JSON STRUCTURE
opsCmd;
param.properties.mat = true;
opsAuth = load(fullfile(gRadar.tmp_path,'ops.mat'));

if opsAuth.isAuthenticated
  param.properties.userName = opsAuth.userName;
  param.properties.isAuthenticated = opsAuth.isAuthenticated;
else
  error('PLEASE LOGIN USING opsLoginUser() FIRST');
end

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
    'Post',{'app' 'rds' 'data' jsonStr 'view' 'getUserProfileData'});
else
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'get/user/profile/data'),'','',...
    'Post',{'app' 'rds' 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,data] = jsonResponseDecode(jsonResponse);

% WRITE OPS AUTH FILE
if status == 1
  message = 'PROFILE WRITTEN SUCCESFULLY';
  save(fullfile(gRadar.tmp_path,'ops.profile.mat'),'-struct','data')
else
  message = data;
end
end