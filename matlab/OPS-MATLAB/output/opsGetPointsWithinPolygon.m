function [ status,data ] = opsGetPointsWithinPolygon( sys, param )
%
% [status,data] = opsGetPointsWithinPolygon(sys,param)
%
% Get the points' info within a polygon
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.location = string ('arctic' or 'antarctic')
%     properties.bound = WKT Polygon Boundary
%     ('POLYGON((-117.56243724291429 -81.19720421323231,-140.87950811907243 -76.39214612404346,-178.98362807884536 -79.05861238803509,173.13100276931294 -86.11254475101619,-117.56243724291429 -81.19720421323231))';)
%     properties.season = string of season list (optional)
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: structure with fields (or error message)
%      properties.Lat = double array
%      properties.Lon = double array
%      properties.Elevation = double array
%      properties.Gps_Time = double array
%      properties.Surface = double array
%      properties.Bottom = double array
%      properties.Thickness = double array
%      properties.Surface_Quality = integer array
%      properties.Bottom_Quality = integer array
%      properties.Season = string cell
%      properties.Frame = string cell
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
    'Post',{'app' sys 'data' jsonStr 'view' 'getPointsWithinPolygon'});
else
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'get/point/polygon'),gOps.dbUser,gOps.dbPswd,...
    'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,decodedJson] = jsonResponseDecode(jsonResponse);

if status == 1
  data.properties.Lat = double(cat(2,decodedJson.Lat{:}));
  data.properties.Lon = double(cat(2,decodedJson.Lon{:}));
  data.properties.Elevation = double(cat(2,decodedJson.Elevation{:}));
  data.properties.Gps_Time = double(cat(2,decodedJson.Gps_Time{:}));
  data.properties.Surface = double(cat(2,decodedJson.Surface{:}));
  data.properties.Bottom = double(cat(2,decodedJson.Bottom{:}));
  data.properties.Thickness = double(cat(2,decodedJson.Thickness{:}));
  data.properties.Surface_Quality = int32(cat(2,decodedJson.Surface_Quality{:}));
  data.properties.Bottom_Quality =  int32(cat(2,decodedJson.Bottom_Quality{:}));
  data.properties.Season = cat(2,decodedJson.Season);
  data.properties.Frame = cat(2,decodedJson.Frame); 
else
  data.properties.Lat = [];
  data.properties.Lon = [];
  data.properties.Elevation = [];
  data.properties.Gps_Time = [];
  data.properties.Surface = [];
  data.properties.Bottom = [];
  data.properties.Thickness = [];
  data.properties.Surface_Quality = [];
  data.properties.Bottom_Quality = [];
  data.properties.Season = [];
  data.properties.Frame = []; 
end

end

