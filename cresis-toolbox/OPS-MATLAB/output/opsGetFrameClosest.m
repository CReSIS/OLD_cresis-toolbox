function [status,data]=opsGetFrameClosest(sys,param)
%
% [status,data] = opsGetFrameClosest(sys,param)
%
% Find the closest frame to a given point from the database.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%       properties.location = string ('arctic' or 'antarctic')
%       properties.x = double
%       properties.y = double
%       OPTIONAL:
%           properties.season = string OR cell of string/s
%           properties.startseg = string of minimum segment name
%           properties.stopseg = string of maximum segment name
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: structure with fields (or error message)
%       properties.frame = string
%       properties.X = double array
%       properties.Y = double array
%       properties.gps_time = double array
%       properties.frame = string
%       properties.season = string
%       properties.segment_id = double
%       properties.start_gps_time = double
%       properties.stop_gps_time = double
%
% Author: Kyle W. Purdon, Trey Stafford

% CONSTRUCT THE JSON STRUCTURE
param.properties.status = [true,false];
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
    'Post',{'app' sys 'data' jsonStr 'view' 'getFrameClosest'});
else
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'get/frame/closest'),gOps.dbUser,gOps.dbPswd,...
    'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,decodedJson] = jsonResponseDecode(jsonResponse);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
data.properties.frame = decodedJson.frame;
data.properties.season = decodedJson.season;
data.properties.segment_id = double(decodedJson.segment_id);
data.properties.start_gps_time = double(decodedJson.start_gps_time);
data.properties.stop_gps_time = double(decodedJson.stop_gps_time);
data.properties.X = cell2mat(decodedJson.X)';
data.properties.Y = cell2mat(decodedJson.Y)';
data.properties.gps_time = cell2mat(decodedJson.gps_time)';
% data.properties.echograms = decodedJson.echograms;

end