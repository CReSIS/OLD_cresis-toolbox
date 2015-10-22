function [status,data] = ops_get_initial_data(sys,param)
%
% [status,data] = ops_get_initial_data(sys,param)
%
% Creates an initial datapack for field depolyments of the OPS system 
% containing all database information associated with segments marked by 
% generic column in one or more param sheets. Datapacks should be placed
% in the ../vagrant/data/postgresql/ directory prior to 'vagrant up.' 
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     param_sheets = string or cell of strings ('..\params-cr1\rds_param_2013_Greenland_P3')
%     OPTIONAL:
%       layers = string or cell of strings. Limits initial data to specified layers ('surface','bottom') 
% 
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: string contianing url to download datapack from.
%   A matlab web-browser prompting the datapack download
%
% Author: Trey Stafford
% 
    %Initialize properties 
    properties.radars={};properties.seasons={};properties.segments={};

    % READ ALL PARAM SHEETS PROVIDED
    for xls_idx = 1:length(param.param_sheets)
        fprintf('READING PARAM SHEET %d OF %d\n',xls_idx,length(param.param_sheets))
        params = read_param_xls(param.param_sheets{xls_idx});
        %GET THE RADAR AND SEASON NAMES FOR THE CURRENT SHEET
        properties.radars(end+1) = {params(1).radar_name};
        properties.seasons(end+1) = {params(1).season_name};
        for param_idx = 1:length(params)
          param_row = params(param_idx);
          if param_row.cmd.generic == 1
            if ~isempty(regexpi(param_row.cmd.notes,'do not process'))
              warning('You have enabled a segment with ''do not process'' in the cmd.notes, dbcont to continue');
              keyboard
            end
            %GET ALL SPECIFIED SEGMENTS 
            properties.segments(end+1) = {param_row.day_seg};
          end
        end
    end
    
    % KEEP ONLY UNIQUE VALUES
    properties.segments = unique(properties.segments);
    properties.radars = unique(properties.radars);
    
    % CHECK IF ONLY SOME LAYERS ARE SPECIFIED
    if isfield(param,'layers')
        properties.layers = param.layers;
    end
    
    % CONSTRUCT THE JSON STRUCTURE
    json_struct = struct('type','Feature','properties',properties);
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
        'Post',{'app' sys 'data' json_str 'view' 'get_initial_data'});
    else
      [json_response,~] = cr_urlread(strcat(server_url,'get/initial/data'),db_user,db_pswd,...
        'Post',{'app' sys 'data' json_str});
    end

% DECODE THE SERVER RESPONSE
[status,decoded_json] = json_response_decode(json_response);

% OPEN A WEB BROWSER TO BEGIN DOWNLOAD.
if status==1
    data = strcat(sys_url,decoded_json);
    web(strcat(sys_url,decoded_json))   
else
    error(decoded_json)

end