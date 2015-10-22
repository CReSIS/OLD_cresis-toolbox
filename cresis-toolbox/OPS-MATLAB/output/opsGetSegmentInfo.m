function [status,data] = opsGetSegmentInfo(sys,param)
%
% [status,data] = opsGetSegmentInfo(sys,param)
%
% Retrieves information for a single segment id from the database.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.segment_id = integer
%         or
%     properties.segment = string
%     properties.season = string
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: structure with fields (or error message)
%       properties.season = string
%       properties.segment = string
%       properties.frame = cell array
%       properties.start_gps_time = double array
%       properties.stop_gps_time = double array
%
% Author: Kyle W. Purdon, Trey Stafford, Weibo Liu

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
        'Post',{'app' sys 'data' jsonStr 'view' 'getSegmentInfo'});
else
    [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'get/segment/info'),gOps.dbUser,gOps.dbPswd,...
        'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,decodedJson] = jsonResponseDecode(jsonResponse);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
data.properties.season = decodedJson.season;
data.properties.segment = decodedJson.segment;
data.properties.frame = decodedJson.frame;
data.properties.start_gps_time = double(cell2mat(decodedJson.start_gps_time)');
data.properties.stop_gps_time = double(cell2mat(decodedJson.stop_gps_time)');

end