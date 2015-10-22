% Creates an initial datapack for field depolyments of the OPS system 
% containing all database information associated with segments marked by 
% generic column in one or more param sheets. Datapacks should be placed
% in the ../vagrant/data/postgresql/ directory prior to 'vagrant up.' 
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param_sheets = string or cell of strings ('..\params-cr1\rds_param_2013_Greenland_P3')
%   OPTIONAL:
%    param.properties.layers = string or cell of strings. Limits initial data to specified layers ('surface','bottom') 
% 
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   data: string contianing url to download datapack from.
%
% See also: opsGetInitialData.m
% 
% Author: Trey Stafford
% 
%% USER INPUT
sys = 'rds';
param_sheets = {'H:\scripts\params-cr1\rds_param_2013_Greenland_P3.xls','H:\scripts\params-cr1\rds_param_2010_Greenland_P3.xls'};
param.properties.layers = {'surface','bottom'};

%% AUTOMATED SECTION
%Initialize properties 
param.properties.radars={};param.properties.seasons={};param.properties.segments={};

% READ ALL PARAM SHEETS PROVIDED
for xls_idx = 1:length(param_sheets)
    fprintf('READING PARAM SHEET %d OF %d\n',xls_idx,length(param_sheets))
    params = read_param_xls(param_sheets{xls_idx});
    %GET THE RADAR AND SEASON NAMES FOR THE CURRENT SHEET
    properties.radars{end+1} = params(1).radar_name;
    param.properties.seasons{end+1} = params(1).season_name;
    for param_idx = 1:length(params)
      param_row = params(param_idx);
      if param_row.cmd.generic == 1
        if ~isempty(regexpi(param_row.cmd.notes,'do not process'))
          warning('You have enabled a segment with ''do not process'' in the cmd.notes, dbcont to continue');
          keyboard
        end
        %GET ALL SPECIFIED SEGMENTS 
        param.properties.segments{end+1} = param_row.day_seg;
      end
    end
end

% KEEP ONLY UNIQUE VALUES
param.properties.segments = unique(param.properties.segments);
param.properties.radars = unique(param.properties.radars);

% RUN opsGetInitialData
[status,data] = opsGetInitialData(sys,param);
if status == 1
    fprintf('Initial datapack successfully created.\n Download at: %s\n', data);
else
    fprintf('Problem Creating Initial Datapack... \n')
    error(data)
end