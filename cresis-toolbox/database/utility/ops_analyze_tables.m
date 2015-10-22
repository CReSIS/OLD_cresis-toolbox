function [status,message] = ops_analyze_tables(sys,param)
%
% [status,message] = ops_analyze_tables(sys,param)
%
% Analyzes tables and updates database statistics.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.tables = cell array of strings specifying table names to be analyzed. (without system_ prefix) 
%
% Output:
%   status: integer (0:Error,1:Success)
%   message: status message
%
% Author: Trey Stafford

% CONSTRUCT THE JSON STRUCTURE
json_struct = struct('type','Feature','properties',param.properties);

% CONVERT THE JSON STRUCUTRE TO A JSON STRING
try
  json_str = tojson(json_struct);
catch ME
  json_str = savejson('',json_struct,'FloatFormat','%2.10f','NaN','null');
end

% SEND THE COMMAND TO THE SERVER
ops_sys_cmd;
if profile_cmd
  [json_response,~] = cr_urlread(strcat(server_url,'profile'),db_user,db_pswd,...
    'Post',{'app' sys 'data' json_str 'view' 'tables_analyze'});
else
  [json_response,~] = cr_urlread(strcat(server_url,'tables/analyze'),db_user,db_pswd,...
    'Post',{'app' sys 'data' json_str});
end

% DECODE THE SERVER RESPONSE
[status,decoded_json] = json_response_decode(json_response);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
message = decoded_json;

end