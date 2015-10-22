function [status,data] = ops_get_closest_point(sys,param)
%
% [status,data] = ops_get_closest_point(sys,param)
%
% Find the closest point to a given point from the database.
%
% Input:
%   sys: (string) system name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.location = string ('arctic' or 'antarctic')
%     properties.x = double
%     properties.y = double
%     properties.start_gps_time = double OR double array
%     properties.stop_gps_time = double OR double array
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: structure with fields (or error message)
%       properties.gps_time = double
%       properties.X = double
%       properties.Y = double
%
% Author: Kyle W. Purdon, Trey Stafford

% CONSTRUCT THE JSON STRUCTURE
json_struct = struct('type','Feature','properties',param.properties);

% CONVERT THE JSON STRUCTURE TO A JSON STRING
try
  json_str = tojson(json_struct);
catch ME
  json_str = savejson('',json_struct,'FloatFormat','%2.10f','NaN','null');
end

% SEND THE COMMAND TO THE SERVER
ops_sys_cmd;
if profile_cmd
  [json_response,~] = cr_urlread(strcat(server_url,'profile'),db_user,db_pswd,...
    'Post',{'app' sys 'data' json_str 'view' 'get_closest_point'});
else
  [json_response,~] = cr_urlread(strcat(server_url,'get/closest/point'),db_user,db_pswd,...
    'Post',{'app' sys 'data' json_str});
end

% DECODE THE SERVER RESPONSE
[status,decoded_json] = json_response_decode(json_response);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
data.properties.X = decoded_json.X;
data.properties.Y = decoded_json.Y;
data.properties.gps_time = decoded_json.gps_time;

end