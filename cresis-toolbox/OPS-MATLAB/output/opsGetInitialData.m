function [status,data_url] = opsGetInitialData(sys,param)
%
% [status,data_url] = opsGetInitialData(sys,param)
%
% Creates an initial datapack for field depolyments of the OPS system
% containing all database information associated with a given set of
% segments. Datapacks should be placed in the
% ../vagrant/data/postgresql/ directory prior to running provisions.sh
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.seasons = cell of strings
%	  properties.segments = cell of strings.
%	  properties.radars = cell of strings.
%     OPTIONAL:
%       properties.layers = string or cell of strings. Limits initial data to specified layers ('surface','bottom')
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data_url: string contianing url to download datapack from.
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
        'Post',{'app' sys 'data' jsonStr 'view' 'getInitialData'});
else
    [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'get/initial/data'),gOps.dbUser,gOps.dbPswd,...
        'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,decodedJson] = jsonResponseDecode(jsonResponse);

if status==1
    data_url = strcat(gOps.sysUrl,decodedJson);
    fprintf('Initial Datapack Created.\n Download from %s \n', data_url);
end