function [status,message] = ops_create_layer_points(sys,param)
%
% [status,message] = ops_create_layer_points(sys,param)
%
% Loads a layer into the database.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     geometry.coordinates = double array of format ([lon; lat; elevation;])
%     properties.username = string
%     properties.location = string ('arctic' or 'antarctic')
%     properties.segment = string
%     properties.gps_time = double array
%     properties.twtt = double array
%     properties.type = integer arry (0,1 or 2)
%     properties.quality = integer array (1,2 or 3)
%     properties.lyr_name = string ('surface','bottom', etc...)
%     properties.lyr_group_name = string (optional)
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   message: status message or data
%
% Author: Kyle W. Purdon

% CONSTRUCT THE JSON STRUCTURE
if numel(param.geometry.coordinates) <= 3
  geom_type = 'Point';
else
  geom_type = 'LineString';
end

% SET UP DEFAULT VALUES FOR GROUP
if ~isfield(param.properties,'lyr_group_name')
  if any(strcmp(param.properties.lyr_name,{'surface','bottom'}))
    param.properties.lyr_group_name = 'standard';
  else
    param.properties.lyr_group_name = 'uncatagorized';
  end
end

% CONSTRUCT THE JSON STRUCTURE
json_struct = struct('type','Feature','geometry',struct('type',geom_type,'coordinates',param.geometry.coordinates'),...
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
    'Post',{'app' sys 'data' json_str 'view' 'create_layer_points'});
else
  [json_response,~] = cr_urlread(strcat(server_url,'create/layer/points'),db_user,db_pswd,...
    'Post',{'app' sys 'data' json_str});
end

% DECODE THE SERVER RESPONSE
[status,decoded_json] = json_response_decode(json_response);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
message = decoded_json;

end