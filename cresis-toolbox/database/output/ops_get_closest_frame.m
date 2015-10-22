function [status,data]=ops_get_closest_frame(sys,param)
%
% [status,data] = ops_get_closest_frame(sys,param)
%
% Find the closest frame to a given point from the database.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.location = string ('arctic' or 'antarctic')
%     properties.season = string OR cell of string/s (can be '')
%     properties.x = double
%     properties.y = double
%     properties.status = string OR cell of string/s (can be '')
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: structure with fields (or error message)
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
param.properties.status = {'public','private'};
json_struct = struct('type','Feature','properties',param.properties);

% CONVERT THE JSON STRUCTURE TO A JSON STRING
try
  json_str = tojson(json_struct);
catch ME
  json_str = savejson('',json_struct,'FloatFormat','%2.10f','NaN','null');
end

% SEND THE COMMAND TO THE SERVER
ops_sys_cmd
if profile_cmd
  [json_response,~] = cr_urlread(strcat(server_url,'profile'),db_user,db_pswd,...
    'Post',{'app' sys 'data' json_str 'view' 'get_closest_frame'});
else
  [json_response,~] = cr_urlread(strcat(server_url,'get/closest/frame'),db_user,db_pswd,...
    'Post',{'app' sys 'data' json_str});
end

% DECODE THE SERVER RESPONSE
[status,decoded_json] = json_response_decode(json_response);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
data.properties.frame = decoded_json.frame;
data.properties.season = decoded_json.season;
data.properties.segment_id = decoded_json.segment_id;
data.properties.start_gps_time = decoded_json.start_gps_time;
data.properties.stop_gps_time = decoded_json.stop_gps_time;
data.properties.X = double(cat(2,decoded_json.X{:}));
data.properties.Y = double(cat(2,decoded_json.Y{:}));
data.properties.gps_time = double(cat(2,decoded_json.gps_time{:}));
data.properties.echogram_url = decoded_json.echogram_url;

% FORCE X,Y TO BE SORTED ON gps_time
[data.properties.gps_time,sort_idxs] = sort(data.properties.gps_time);
data.properties.X = data.properties.X(sort_idxs);
data.properties.Y = data.properties.Y(sort_idxs);

end