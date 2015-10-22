function [status,data] = ops_get_path(sys,param)
%
% [status,data] = ops_get_path(sys,param)
%
% Gets the point path from the database based on the given params
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.location = string ('arctic' or 'antarctic')
%     properties.season = string
%     properties.start_gps_time = double
%     properties.stop_gps_time = double
%
% Output:
%   status: integer (0:Error,1:Success)
%   data: structure with fields (or error message)
%       properties.gps_time = double array
%       properties.X = double array
%       properties.Y = double array
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
    'Post',{'app' sys 'data' json_str 'view' 'get_path'});
else
  [json_response,~] = cr_urlread(strcat(server_url,'get/path'),db_user,db_pswd,...
    'Post',{'app' sys 'data' json_str});
end

% DECODE THE SERVER RESPONSE
[status,decoded_json] = json_response_decode(json_response);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
data.properties.gps_time = double(cat(2,decoded_json.gps_time{:}));
data.properties.X = double(cat(2,decoded_json.X{:}));
data.properties.Y = double(cat(2,decoded_json.Y{:}));

end