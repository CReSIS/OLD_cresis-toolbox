function [status,data] = opsAlterPathResolution(sys,param)
%
% [status,message] = opsAlterPathResolution(sys,param)
%
% Alters a segment's resolution (vertex point spacing along line).
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.segment: {string(s)} OR segment_id [integer(s)]
%     properties.season: (string) the season name of the segment(s) to be altered.
%     properties.resolution: number of meters between each point
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: status message
%
% Author: Trey Stafford

% CONSTRUCT THE JSON STRUCTURE
[param,~,~] = opsAuthenticate(param);
jsonStruct = struct('properties',param.properties);

% CONVERT THE JSON STRUCUTRE TO A JSON STRING
try
  jsonStr = tojson(jsonStruct);
catch ME
  jsonStr = savejson('',jsonStruct,'FloatFormat','%2.10f','NaN','null');
end

% SEND THE COMMAND TO THE SERVER
opsCmd;
if gOps.profileCmd
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'profile'),gOps.dbUser,gOps.dbPswd,...
    'Post',{'app' sys 'data' jsonStr 'view' 'alterPathResolution'});
else
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'alter/path/resolution'),gOps.dbUser,gOps.dbPswd,...
    'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,data] = jsonResponseDecode(jsonResponse);

end

