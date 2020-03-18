function [status,data] = opsGetCrossovers(sys,param)
%
% [status,data] = opsGetCrossovers(sys,param)
%
% Retrieves crossovers from the database.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.location = string ('arctic' or 'antarctic')
%     properties.lyr_name = string or cell of strings ('surface','bottom')
%     PICK ONE OF THE FOLLOWING:
%       properties.point_path_id: integer or array of integers
%       OR
%       properties.frame = string or cell of strings
%       properties.segment_id = integer or cell of integers
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: structure with fields (or error message)
%       properties.source_point_path_id = integer array
%       properties.cross_point_path_id = integer array
%       properties.source_elev = double array
%       properties.cross_elev = double array
%       properties.layer_id = integer array
%       properties.frame_id = integer array
%       properties.segment_id = integer array
%       properties.twtt = double array
%       properties.angle = double array
%       properties.abs_error = double array
%       properties.cross_quality = integer array
%
% Author: Trey Stafford, Kyle W. Purdon

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
    'Post',{'app' sys 'data' jsonStr 'view' 'getCrossovers'});
else
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'get/crossovers'),gOps.dbUser,gOps.dbPswd,...
    'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,decodedJson] = jsonResponseDecode(jsonResponse);

if status == 2
  data.properties.source_point_path_id = [];
  data.properties.cross_point_path_id = [];
  data.properties.source_elev = [];
  data.properties.cross_elev = [];
  data.properties.layer_id = [];
  data.properties.frm_str = {};
  data.properties.season_name = {};
  data.properties.seg_id = [];
  data.properties.twtt = [];
  data.properties.angle = [];
  data.properties.abs_error = [];
  data.properties.cross_quality = [];
else
  data.properties.source_point_path_id = cell2mat(decodedJson.source_point_path_id);
  data.properties.cross_point_path_id = cell2mat(decodedJson.cross_point_path_id);
  data.properties.source_elev = cell2mat(decodedJson.source_elev);
  data.properties.cross_elev = cell2mat(decodedJson.cross_elev);
  data.properties.layer_id = cell2mat(decodedJson.layer_id);
  data.properties.frm_str = decodedJson.frame_name;
  data.properties.season_name = decodedJson.season_name;
  data.properties.seg_id = cell2mat(decodedJson.segment_id);
  data.properties.twtt = cell2mat(decodedJson.twtt);
  data.properties.angle = cell2mat(decodedJson.angle);
  data.properties.abs_error = cell2mat(decodedJson.abs_error);
  data.properties.cross_quality = cellfun(@double,decodedJson.cross_quality);
end
end