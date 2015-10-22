function [status,data] = ops_get_layers(sys)
%
% [status,data] = ops_get_layers(sys)
%
% Retrieves all layers for a given system with normal status from the database.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: structure with fields (or error message)
%       properties.lyr_id = double array
%       properties.lyr_name = cell array
%
% Author: Kyle W. Purdon, Trey Stafford

% SEND THE COMMAND TO THE SERVER
ops_sys_cmd;
if profile_cmd
  [json_response,~] = cr_urlread(strcat(server_url,'profile'),db_user,db_pswd,...
    'Post',{'app' sys 'data' '{"nodata":"novalue"}' 'view' 'get_layers'});
else
  [json_response,~] = cr_urlread(strcat(server_url,'get/layers'),db_user,db_pswd,...
    'Post',{'app' sys 'data' '{"nodata":"novalue"}'});
end

% DECODE THE SERVER RESPONSE
[status,decoded_json] = json_response_decode(json_response);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
data.properties.lyr_name = decoded_json.lyr_name;
data.properties.lyr_id = double(cat(2,decoded_json.lyr_id{:}));

end