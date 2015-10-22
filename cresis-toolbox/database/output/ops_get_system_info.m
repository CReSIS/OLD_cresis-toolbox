function [status,data] = ops_get_system_info()
%
% [status,data] = ops_get_system_info()
%
% Retrieves all systems/seasons/locations from the database
%
% Input:
%   none
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: structure with fields (or error message)
%       properties.systems = cell of string/s
%       properties.seasons = cell of string/s
%       properties.locations = cell of string/s
%
% Author: Kyle W. Purdon

% SEND THE COMMAND TO THE SERVER
ops_sys_cmd;
if profile_cmd
  [json_response,~] = cr_urlread(strcat(server_url,'profile'),db_user,db_pswd,...
    'Post',{'view' 'get_system_info'});
else
  [json_response,~] = cr_urlread(strcat(server_url,'get/system/info'),db_user,db_pswd,...
    'Post',{});
end

% DECODE THE SERVER RESPONSE
[status,decoded_json] = json_response_decode(json_response);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
data.properties.systems = decoded_json.systems;
data.properties.seasons = decoded_json.seasons;
data.properties.locations = decoded_json.locations;

end