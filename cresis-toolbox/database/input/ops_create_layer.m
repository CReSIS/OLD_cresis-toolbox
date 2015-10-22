function [status,data] = ops_create_layer(sys,param)
%
% [status,data] = ops_create_layer(sys,param)
%
% Adds a new layer to the database. Restores status if it was deleted.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.lyr_name = string
%     properties.lyr_group_name = string (optional)
%     properties.lyr_description = string (optional)
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: structure with fields (or error/warning message)
%       properties.lyr_id = double
%
% Author: Kyle W. Purdon

% SET UP DEFAULT VALUES FOR GROUP AND DESCRIPTION
if ~isfield(param.properties,'lyr_group_name')
  if any(strcmp(param.properties.lyr_name,{'surface','bottom'}))
    param.properties.lyr_group_name = 'standard';
  else
    param.properties.lyr_group_name = 'uncatagorized';
  end
end
if ~isfield(param.properties,'lyr_description')
  param.properties.lyr_description = '';
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
    'Post',{'app' sys 'data' json_str 'view' 'create_layer'});
else
  [json_response,~] = cr_urlread(strcat(server_url,'create/layer'),db_user,db_pswd,...
    'Post',{'app' sys 'data' json_str});
end

% DECODE THE SERVER RESPONSE
[status,decoded_json] = json_response_decode(json_response);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
data.properties.lyr_id = decoded_json.lyr_id;
data.properties.lyr_group_id = decoded_json.lyr_group_id;

end