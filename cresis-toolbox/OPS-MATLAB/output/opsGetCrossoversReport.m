function [status,data] = opsGetCrossoversReport(sys,param)
%
% [status,data] = opsGetCrossoversReport(sys,param)
%
% Creates a crossovers report (.csv) from crossovers in the database.
% Currently returns only crossovers with reportable errors. 
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.location = string ('arctic' or 'antarctic')
%     properties.lyr_name = string or cell of strings ('surface','bottom')
%       'all' will fetch crossovers for all layers
%     PICK ONE OF THE FOLLOWING:
%       properties.seasons: string or cell of strings (The seasons for
%           which crossovers will be returned. Specifiying one season 
%            will return self-intersecting crossovers for that season only.)
%       properties.frame = string or cell of strings
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: url to crossover report .csv on the server
%
% Author: Trey Stafford, Kyle W. Purdon

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
        'Post',{'app' sys 'data' jsonStr 'view' 'getCrossoversReport'});
else
    [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'get/crossovers/report'),gOps.dbUser,gOps.dbPswd,...
        'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,decodedJson] = jsonResponseDecode(jsonResponse);

if status == 2
    data = []; % CLIENT NEEDS TO HANDLE EMPTY DATA
else
    data = strcat(gOps.sysUrl,decodedJson);
end
end