function [status,data] = ops_bulk_delete(sys,param)
%
% [status,data] = ops_bulk_delete(sys,param)
%
% Completely removes data from OPS. 
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.type = string (type of deletion to perform)
%       'full': deletes all layerdata and path information
%       'layer': deletes only layerdata
%       note: If 'full' and all segments for the season are deleted,
%          the season information will also be removed from the database.
%      properties.segments = cell of segment name(s)
%       examples:
%           {} : deletes all segments in the database for the given season
%           {'20110331_01'}: deletes the given single segment
%           {'20110331_01','20110331_02'}: deletes both given segments
%      properties.season_name: (string) season name for given segments ('2011_Greenland_P3')
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: Message indicating deletion status
%
% Example:
%   [status,data] = ops_bulk_delete('rds',struct('properties',struct('season','2009_Greenland_TO','type','full')));
%   [status,data] = ops_bulk_delete('rds',struct('properties',struct('season','2009_Greenland_TO','type','full','segments',{{'20090402_01','20090402_02'}})));
%
% Author: Trey Stafford, Kyle Purdon 
%      
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
        'Post',{'app' sys 'data' json_str 'view' 'bulk_delete'});
    else
      [json_response,~] = cr_urlread(strcat(server_url,'delete/bulk'),db_user,db_pswd,...
        'Post',{'app' sys 'data' json_str});
    end

% DECODE THE SERVER RESPONSE
[status,decoded_json] = json_response_decode(json_response);

% Return the status of the deletion. 
if status==1
    data = decoded_json;
    fprintf('%s\n', decoded_json);
else
    data = decoded_json;
    error(decoded_json);
end