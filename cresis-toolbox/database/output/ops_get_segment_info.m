function [status,data] = ops_get_segment_info(sys,param)
%
% [status,data] = ops_get_segment_info(sys,param)
%
% Retrieves information for a single segment id from the database.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.location = string ('arctic' or 'antarctic')
%     properties.segment_id = integer
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: structure with fields (or error message)
%       properties.season = string
%       properties.segment = string
%       properties.frame = cell array
%       properties.start_gps_time = double array
%       properties.stop_gps_time = double array
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
    'Post',{'app' sys 'data' json_str 'view' 'get_segment_info'});
else
  [json_response,~] = cr_urlread(strcat(server_url,'get/segment/info'),db_user,db_pswd,...
    'Post',{'app' sys 'data' json_str});
end

% DECODE THE SERVER RESPONSE
[status,decoded_json] = json_response_decode(json_response);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
data.properties.season = decoded_json.season;
data.properties.segment = decoded_json.segment;
data.properties.frame = decoded_json.frame;
data.properties.start_gps_time = double(cat(2,decoded_json.start_gps_time{:}));
data.properties.stop_gps_time = double(cat(2,decoded_json.stop_gps_time{:}));

end