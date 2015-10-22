function [status,data] = ops_get_crossover_error(sys,param)
%
% [status,data] = ops_get_crossover_error(sys,param)
%
% Retrieves crossovers from the database.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.location = string ('arctic' or 'antarctic')
%     properties.lyr_name = string or cell of strings ('surface','bottom')
%     properties.search_offset = (optional) default = 5 seconds
%     PICK ONE OF THE FOLLOWING:
%       properties.season = string or cell of strings
%       properties.segment = string or cell of strings
%       properties.frame = string or cell of strings
%     OPTIONAL:
%       properties.excluded_seasons = string or cell of strings. Excludes seasons from output. 
% 
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: structure with fields (or error message)
%       properties.X = double array
%       properties.Y = double array
%       properties.twtt_1 = double array
%       properties.twtt_2 = double array
%       properties.abs_error = double array
%       properties.cross_angle = double array
%       properties.gps_time_1 = double array
%       properties.gps_time_2 = double array
%       properties.ELEV_1 = double array
%       properties.ELEV_2 = double array
%       properties.distance_1 = double array
%       properties.distance_2 = double array
%       properties.frame_id_1 = integer array
%       properties.frame_id_2 = integer array
%       properties.frame_1 = cell of string/s
%       properties.frame_2 = cell of string/s
%       properties.lyr_id = integer array
%
% Author: Trey Stafford, Kyle W. Purdon

% SET THE DEFAULT SEARCH RADIUS
if ~isfield(param.properties,'search_offset')
  param.properties.search_offset = 5;
end

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
    'Post',{'app' sys 'data' json_str 'view' 'get_crossover_error'});
else
  [json_response,~] = cr_urlread(strcat(server_url,'get/crossover/error'),db_user,db_pswd,...
    'Post',{'app' sys 'data' json_str});
end

% DECODE THE SERVER RESPONSE
[status,decoded_json] = json_response_decode(json_response);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
if status == 2
  data.properties.X = [];
  data.properties.Y = [];
  data.properties.twtt_1 = [];
  data.properties.twtt_2 = [];
  data.properties.abs_error = [];
  data.properties.cross_angle = [];
  data.properties.gps_time_1 = [];
  data.properties.gps_time_2 = [];
  data.properties.ELEV_1 = [];
  data.properties.ELEV_2 = [];
  data.properties.distance_1 = [];
  data.properties.distance_2 = [];
  data.properties.frame_id_1 = [];
  data.properties.frame_id_2 = [];
  data.properties.frame_1 = {};
  data.properties.frame_2 = {};
  data.properties.lyr_id = [];
else
  data.properties.X = double(cat(2,decoded_json.X{:}));
  data.properties.Y = double(cat(2,decoded_json.Y{:}));
  data.properties.twtt_1 = double(cat(2,decoded_json.twtt_1{:}));
  data.properties.twtt_2 = double(cat(2,decoded_json.twtt_2{:}));
  data.properties.abs_error = double(cat(2,decoded_json.abs_error{:}));
  data.properties.cross_angle = double(cat(2,decoded_json.cross_angle{:}));
  data.properties.gps_time_1 = double(cat(2,decoded_json.gps_time_1{:}));
  data.properties.gps_time_2 = double(cat(2,decoded_json.gps_time_2{:}));
  data.properties.ELEV_1 = double(cat(2,decoded_json.ELEV_1{:}));
  data.properties.ELEV_2 = double(cat(2,decoded_json.ELEV_2{:}));
  data.properties.distance_1 = double(cat(2,decoded_json.distance_1{:}));
  data.properties.distance_2 = double(cat(2,decoded_json.distance_2{:}));
  data.properties.frame_id_1 = double(cat(2,decoded_json.frame_id_1{:}));
  data.properties.frame_id_2 = double(cat(2,decoded_json.frame_id_2{:}));
  data.properties.frame_1 = decoded_json.frame_1;
  data.properties.frame_2 = decoded_json.frame_2;
  data.properties.lyr_id = double(cat(2,decoded_json.lyr_id{:}));
end
end