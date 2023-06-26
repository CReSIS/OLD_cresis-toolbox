function [status,data] = opsGetFrameSearch(sys,param)
%
% [status,data] = opsGetFrameSearch(sys,param)
%
% Retrieves the closest frame based on a search string from the database.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.search_str = string (e.g., '2010','201005','20100510_04_010',...)
%     properties.location = string ('arctic' or 'antarctic')
%
%     optional:
%       properties.season = string (eg. '2011_Greenland_P3')
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: structure with fields (or error message)
%       properties.season = string
%       properties.segment_id = int
%       properties.frame = string
%       properties.X = double (x coordinate in map projection)
%       properties.Y = double (y coordinate in map projection)
%       properties.gps_time = double (z coordinate)
%

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
    'Post',{'app' sys 'data' jsonStr 'view' 'getFrameSearch'});
else
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'get/frame/search'),gOps.dbUser,gOps.dbPswd,...
    'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,decodedJson] = jsonResponseDecode(jsonResponse);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
if status == 2
  data = [];
else
  data.properties.frame = decodedJson.frame;
  data.properties.season = decodedJson.season;
  data.properties.segment_id = double(decodedJson.segment_id);
  data.properties.X = cell2mat(decodedJson.X)';
  data.properties.Y = cell2mat(decodedJson.Y)';
  data.properties.gps_time = cell2mat(decodedJson.gps_time)';
end
end