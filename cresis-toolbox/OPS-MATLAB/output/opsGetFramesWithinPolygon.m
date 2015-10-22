function [ status,message ] = opsGetFramesWithinPolygon( sys, param )
%
% [status,message] = opsGetFramesWithinPolygon(sys,param)
%
% Get the frame names within a polygon
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.location = string ('arctic' or 'antarctic')
%     properties.bound = WKT Polygon Boundary
%     ('POLYGON((-117.56243724291429 -81.19720421323231,-140.87950811907243 -76.39214612404346,-178.98362807884536 -79.05861238803509,173.13100276931294 -86.11254475101619,-117.56243724291429 -81.19720421323231))';)
%     properties.season = string of season list (optional)
% Output:
%   status: integer (0:Error,1:Success, 2:Warning)
%   message: status message
%
% Author: Weibo Liu

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
    'Post',{'app' sys 'data' jsonStr 'view' 'get/frame/polygon'});
else
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'get/frame/polygon'),gOps.dbUser,gOps.dbPswd,...
    'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,decodedJson] = jsonResponseDecode(jsonResponse);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
message = decodedJson;

end

