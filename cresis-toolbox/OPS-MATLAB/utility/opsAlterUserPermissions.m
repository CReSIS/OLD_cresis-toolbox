function [status,data] = opsAlterUserPermissions(sys,param)
%
% [status,data] = opsAlterUserPermissions(sys,param)
%
% Alters a user's OPS database permissions.
%
% Input: 
%   sys: (string) sys name ('rds','accum','snow',...)
%       (sys only used in altering layer and/or season groups)
%   param: structure with fields
%       user_name: (string) the username of the user for which permissions
%           will be altered
%       One or more of the following (cell array with '+' appends permissions):
%           sys_layer_groups: integer(s) or cell array of integer(s) with '+'
%           sys_season_groups: integer(s) or cell array of integer(s) with '+'
%           layerGroupRelease: (string) 'True' or 'False'
%           bulkDeleteData: (string) 'True' or 'False'
%           createData: (string) 'True' or 'False'
%           seasonRelease: (string) 'True' or 'False'
%           isRoot: (string) 'True' or 'False'
%       
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data:(string) A status message 
%
% Author: Trey Stafford
[param,~,opsProfile] = opsAuthenticate(param,true);
if ~opsProfile.isRoot
    error('USER NOT AUTHORIZED TO ALTER PERMISSIONS');
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
opsCmd;
if gOps.profileCmd
    [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'profile'),gOps.dbUser,gOps.dbPswd,...
        'Post',{'app' sys 'data' jsonStr 'view' 'alterUserPermissions'});
else
    [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'alter/user/permissions'),gOps.dbUser,gOps.dbPswd,...
        'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,data] = jsonResponseDecode(jsonResponse);

if status == 1
    fprintf('%s \n',data)
end

end

