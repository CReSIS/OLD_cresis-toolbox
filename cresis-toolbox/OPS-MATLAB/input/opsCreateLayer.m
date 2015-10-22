function [status,data] = opsCreateLayer(sys,param)
%
% [status,data] = opsCreateLayer(sys,param)
%
% Adds a new layer to the database. Restores status if it was deleted.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.lyr_name = string
%     properties.lyr_group_name = string (optional)
%     properties.lyr_description = string (optional)
%     properties.public = True (optional)
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: structure with fields (or error/warning message)
%       properties.lyr_id = integer
%       properties.lyr_group_id = integer
%
% Author: Kyle W. Purdon

% SET UP DEFAULT VALUES FOR GROUP AND DESCRIPTION
if ~isfield(param.properties,'lyr_group_name') || isempty(param.properties.lyr_group_name)
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
[param,~,~] = opsAuthenticate(param);
jsonStruct = struct('properties',param.properties);

% CONVERT THE JSON STRUCTURE TO A JSON STRING
try
  jsonStr = tojson(jsonStruct);
catch ME
  jsonStr = savejson('',jsonStruct,'FloatFormat','%2.10f','NaN','null');
end

% SEND THE COMMAND TO THE SERVER
opsCmd;
if gOps.profileCmd
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'profile'),gOps.dbUser,gOps.dbPswd,...
    'Post',{'app' sys 'data' jsonStr 'view' 'createLayer'});
else
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'create/layer'),gOps.dbUser,gOps.dbPswd,...
    'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,decoded_json] = jsonResponseDecode(jsonResponse);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
data.properties.lyr_id = decoded_json.lyr_id;
data.properties.lyr_group_id = decoded_json.lyr_group_id;

end