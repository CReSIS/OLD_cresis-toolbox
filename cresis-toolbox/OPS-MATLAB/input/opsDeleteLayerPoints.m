function [status,message] = opsDeleteLayerPoints(sys,param)
%
% [status,message] = opsDeleteLayerPoints(sys,param)
%
% Removes layer points from the database.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%
%     properties.start_point_path_id = integer
%     properties.stop_point_path_id = integer
%     properties.max_twtt = double
%     properties.min_twtt = double
%     properties.lyr_name OR properties.lyr_id
%       lyr_name: string ('surface','bottom', etc...)
%       lyr_id = scalar integer (database ID)
%
%     OR
%
%     properties.start_gps_time = double
%     properties.stop_gps_time = double
%     properties.max_twtt = double
%     properties.min_twtt = double
%     properties.lyr_name OR properties.lyr_id
%       lyr_name: string ('surface','bottom', etc...)
%       lyr_id = scalar integer (database ID)
%     properties.segment = string
%     properties.season = string
%     properties.location = string ('arctic' or 'antarctic')
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   message: status message
%
% Author: Kyle W. Purdon, Weibo Liu

% CONSTRUCT THE JSON STRUCTURE
jsonStruct = struct('properties',param.properties);

% CONVERT THE JSON STRUCUTRE TO A JSON STRING
try
    jsonStr = tojson(jsonStruct);
catch ME
    jsonStr = savejson('',jsonStruct,'FloatFormat','%2.10f','NaN','null');
end

% SEND THE COMMAND TO THE SERVER
opsCmd;
if gOps.profileCmd
    [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'profile'),gOps.dbUser,gOps.dbPswd,...
        'Post',{'app' sys 'data' jsonStr 'view' 'deleteLayerPoints'});
else
    [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'delete/layer/points'),gOps.dbUser,gOps.dbPswd,...
        'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,message] = jsonResponseDecode(jsonResponse);

end