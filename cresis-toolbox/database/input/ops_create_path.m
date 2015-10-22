function [status,message] = ops_create_path(sys,param)
%
% [status,message] = ops_create_path(sys,param)
%
% Creates a path in the database.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     geometry.coordinates = double array of format ([lon lat gps_time])
%     properties.location = string ('arctic' or 'antarctic')
%     properties.season = string
%     properties.radar = string
%     properties.segment = string
%     properties.segment_start_gps_time = double
%     properties.segment_stop_gps_time = double
%     properties.elev = double array
%     properties.roll = double array
%     properties.pitch = double array
%     properties.heading = double array
%     properties.frame_count = integer array
%     properties.frame_start_gps_time = double array
%     properties.frame_stop_gps_time = double array
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   message: status message or data
%
% Author: Kyle W. Purdon

% CONSTRUCT THE JSON STRUCTURE
json_struct = struct('type','Feature','geometry',struct('type','LineString','coordinates',param.geometry.coordinates'),...
  'properties',param.properties);

% CONVERT THE JSON STRUCUTRE TO A JSON STRING
try
  json_str = tojson(json_struct);
catch ME
  json_str = savejson('',json_struct,'FloatFormat','%2.10f','NaN','null');
end

% SEND THE COMMAND TO THE SERVER
ops_sys_cmd;
if profile_cmd
  [json_response,~] = cr_urlread(strcat(server_url,'profile'),db_user,db_pswd,...
    'Post',{'app' sys 'data' json_str 'view' 'create_path'});
else
  [json_response,~] = cr_urlread(strcat(server_url,'create/path'),db_user,db_pswd,...
    'Post',{'app' sys 'data' json_str});
end

% DECODE THE SERVER RESPONSE
 [status,decoded_json] = json_response_decode(json_response);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
message = decoded_json;

end