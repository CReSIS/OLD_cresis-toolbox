function [status,data] = ops_get_layer_points(sys,param)
%
% [status,data] = ops_get_layer_points(sys,param)
%
% Gets the layer points from the database for the given params.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.location = string ('arctic' or 'antarctic')
%     properties.season = string OR cell of string/s
%
%     select one:
%       properties.segment_id = integer
%       properties.segment = string
%
%     optional (will use entire segment gps_time range if not given)
%       properties.start_gps_time = double
%       properties.stop_gps_time = double
%       properties.return_geom = string ('none','geog','proj')
%         none returns no geometry (default)
%         geog returns lat/lon/elev in wsg84 (degrees)
%         proj returns x/y/elev in polar stereographic (meters)
%
%     properties.lyr_name = string or cell of strings ('surface' OR {'suface','bottom'} OR 'all')
%       'all' : all layers present in the database
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: structure with fields (or error message)
%       properties.lyr_id = integer array
%       properties.gps_time = double array
%       properties.twtt = double array
%       properties.type = integer array
%       properties.quality = integer array
%
%       optional:
%         properties.
%
% Author: Kyle W. Purdon

% SET UP DEFAULTS
if ~isfield(param.properties,'return_geom') 
  param.properties.return_geom = 'none';
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
    'Post',{'app' sys 'data' json_str 'view' 'get_layer_points'});
else
  [json_response,~] = cr_urlread(strcat(server_url,'get/layer/points'),db_user,db_pswd,...
    'Post',{'app' sys 'data' json_str});
end

% DECODE THE SERVER RESPONSE
[status,decoded_json] = json_response_decode(json_response);

% CREATE THE DATA OUTPUT STRUCTURE OR MESSAGE
if status == 2
  % CREATE EMPTY VARIABLES IF WARNING
  data.properties.gps_time = [];
  data.properties.twtt = [];
  data.properties.type = [];
  data.properties.quality = [];
  data.properties.lyr_id = [];
  if strcmpi(param.properties.return_geom,'geog')
    data.properties.lat = [];
    data.properties.lon = [];
    data.properties.elev = [];
  elseif strcmpi(param.properties.return_geom,'proj')
    data.properties.x = [];
    data.properties.y = [];
    data.properties.elev = [];
  end
else
  % STORE THE UNIQUE VARIABLES IN DATA
  data.properties.gps_time = double(cat(2,decoded_json.gps_time{:}));
  data.properties.twtt = double(cat(2,decoded_json.twtt{:}));
  data.properties.type = double(cat(2,decoded_json.type{:}));
  data.properties.quality = double(cat(2,decoded_json.quality{:}));
  data.properties.lyr_id = double(cat(2,decoded_json.lyr_id{:}));
  if strcmpi(param.properties.return_geom,'geog')
    data.properties.lat = double(cat(2,decoded_json.lat{:}));
    data.properties.lon = double(cat(2,decoded_json.lon{:}));
    data.properties.elev = double(cat(2,decoded_json.elev{:}));
  elseif strcmpi(param.properties.return_geom,'proj')
    data.properties.x = [];
    data.properties.y = [];
    data.properties.elev = [];
  end
end
end