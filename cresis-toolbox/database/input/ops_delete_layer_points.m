function [status,message] = ops_delete_layer_points(sys,param)
%
% [status,message] = ops_delete_layer_points(sys,param)
%
% Removes layer points from the database.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.location = string ('arctic' or 'antarctic')
%     properties.segment = string
%     properties.start_gps_time = double
%     properties.stop_gps_time = double
%     properties.max_twtt = double
%     properties.min_twtt = double
%     properties.lyr_name = string
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   message: status message or data
%
% Author: Kyle W. Purdon

% CONSTRUCT THE JSON STRUCTURE
json_struct = struct('type','Feature','properties',param.properties);

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
    'Post',{'app' sys 'data' json_str 'view' 'delete_layer_points'});
else
  [json_response,~] = cr_urlread(strcat(server_url,'delete/layer/points'),db_user,db_pswd,...
    'Post',{'app' sys 'data' json_str});
end

% DECODE THE SERVER RESPONSE
[status,decoded_json] = json_response_decode(json_response);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
message = decoded_json;

end