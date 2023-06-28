function [status,data] = opsDeleteLayer(sys,param)
%
% [status,data] = opsDeleteLayer(sys,param)
%
% Sets the status of a layer to deleted.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.lyr_name = string
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: structure with fields (or error/warning message)
%       properties.lyr_id = integer
%
% Author: Kyle W. Purdon

% CONSTRUCT THE JSON STRUCTURE
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
        'Post',{'app' sys 'data' jsonStr 'view' 'deleteLayer'});
else
    [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'delete/layer'),gOps.dbUser,gOps.dbPswd,...
        'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,decoded_json] = jsonResponseDecode(jsonResponse);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
data.properties.lyr_id = decoded_json.lyr_id;

end