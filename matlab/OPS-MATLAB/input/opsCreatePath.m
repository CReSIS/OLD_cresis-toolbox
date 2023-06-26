function [status,message] = opsCreatePath(sys,param)
%
% [status,message] = opsCreatePath(sys,param)
%
% Creates a path in the database.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     geometry.coordinates = double array of format ([lon lat])
%     properties.location = string ('arctic' or 'antarctic')
%     properties.season = string
%     properties.radar = string
%     properties.segment = string
%     properties.gps_time = double array
%     properties.elev = double array
%     properties.roll = double array
%     properties.pitch = double array
%     properties.heading = double array
%     properties.frame_count = integer array
%     properties.frame_start_gps_time = double array
%
% Optional Input:
%     properties.season_group = string (default = 'cresis_private')
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   message: status message or data
%
% Author: Kyle W. Purdon

% SET SEASON_GROUP DEFAULTS
if ~isfield(param.properties,'season_group')
  param.properties.season_group = 'cresis_private';
end

% CONSTRUCT THE JSON STRUCTURE
[param,~,~] = opsAuthenticate(param);
jsonStruct = struct('type','Feature','geometry',struct('type','LineString','coordinates',param.geometry.coordinates'),...
    'properties',param.properties);

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
        'Post',{'app' sys 'data' jsonStr 'view' 'createPath'});
else
    [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'create/path'),gOps.dbUser,gOps.dbPswd,...
        'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,message] = jsonResponseDecode(jsonResponse);

end