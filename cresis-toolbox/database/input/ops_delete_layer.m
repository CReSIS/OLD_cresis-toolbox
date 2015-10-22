function [status,data] = ops_delete_layer(sys,param)
%
% [status,data] = ops_delete_layer(sys,param)
%
% Sets the status of a layer to deleted.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.location = string ('arctic' or 'antarctic')
%     properties.lyr_name = string
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: structure with fields (or error/warning message)
%       properties.lyr_id = double
%
% Author: Kyle W. Purdon

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
    'Post',{'app' sys 'data' json_str 'view' 'delete_layer'});
else
  [json_response,~] = cr_urlread(strcat(server_url,'delete/layer'),db_user,db_pswd,...
    'Post',{'app' sys 'data' json_str});
end

% DECODE THE SERVER RESPONSE
[status,decoded_json] = json_response_decode(json_response);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
data.properties.lyr_id = decoded_json.lyr_id;

end