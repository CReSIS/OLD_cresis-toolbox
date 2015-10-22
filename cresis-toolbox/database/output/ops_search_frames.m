function [status,data] = ops_search_frames(sys,param)
%
% [status,data] = ops_search_frames(sys,param)
%
% Retrieves the closest frame based on a search string from the database.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.search_str = string (e.g., '20100510_04_010')
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
%       properties.start_gps_time = double
%       properties.stop_gps_time = double
%       properties.frame = string
%       properties.X = double (x coordinate in map projection)
%       properties.Y = double (y coordinate in map projection)
%       properties.gps_time = double (z coordinate in map projection)
%

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
    'Post',{'app' sys 'data' json_str 'view' 'search_frames'});
else
  [json_response,~] = cr_urlread(strcat(server_url,'search/frames'),db_user,db_pswd,...
    'Post',{'app' sys 'data' json_str});
end

% DECODE THE SERVER RESPONSE
[status,decoded_json] = json_response_decode(json_response);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
if status == 2
  data = [];
else
  data.properties.frame = decoded_json.frame;
  data.properties.season = decoded_json.season;
  data.properties.segment_id = decoded_json.segment_id;
  data.properties.start_gps_time = decoded_json.start_gps_time;
  data.properties.stop_gps_time = decoded_json.stop_gps_time;
  data.properties.X = double(cat(2,decoded_json.X{:}));
  data.properties.Y = double(cat(2,decoded_json.Y{:}));
  data.properties.gps_time = double(cat(2,decoded_json.gps_time{:}));
  
  % FORCE X,Y TO BE SORTED ON gps_time
  [data.properties.gps_time,sort_idxs] = sort(data.properties.gps_time);
  data.properties.X = data.properties.X(sort_idxs);
  data.properties.Y = data.properties.Y(sort_idxs);

end
end