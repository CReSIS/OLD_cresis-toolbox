
%=========================================================================
%
% Using the generic field in the command worksheet of the param spreadsheet, check if segments/frames have been inserted and specified layers are picked. Report on status of each.
%
%  To use this script, manually edit the input fileds under USER INPUT.
%
% Input:
%   sysName: (string) sys name ('rds','accum','snow',...)
%   paramFn: FILENAME.xls OF EXCEL CReSIS PARAMS SHEET
%   location = string ('arctic' or 'antarctic')
%    
%
% Output:
%   report: segments/frames have been inserted or not/ specified layers are picked or not
%
% Author: Weibo Liu
%
% =========================================================================

%% USER INPUT

% Initialize settings structure
settings = [];

% paramFn: FILENAME.xls OF EXCEL CReSIS PARAMS SHEET
settings.paramFn = 'rds_param_2013_Antarctica_P3.xls';

% ----------------------------------------------------------------
% location: LOCATION NAME ('arctic' OR 'antarctic')
% settings.location = 'arctic';
settings.location = 'antarctic';

% ----------------------------------------------------------------
% sysName: SYSTEM NAME ('rds','snow','accum','kuband')
settings.sysName = 'rds';


%% AUTOMATED SECTION
% ERROR CHECK THE INPUT
if ~any(strcmp(settings.location,{'arctic','antarctic'}));
  tmp = questdlg({'Location must be ''arctic'' or ''antarctic''','',sprintf('You entered: %s',settings.location),''},'INVALID LOCATION','OK','OK');
  error('INVALID LOCATION');
elseif ~any(strcmp(settings.sysName,{'rds','snow','accum','kuband'}));
  tmp = questdlg({'System must be ''rds'',''snow'',''accum'' or ''kuband''','',sprintf('You entered: %s',settings.sysName),''},'INVALID SYSTEM','OK','OK');
  error('INVALID SYSTEM');
end

% GET THE CReSIS GLOBAL
global gRadar;
if isvarname('gRadar')
  settings = mergestruct(settings,gRadar);
else
  warning('gRadar IS NOT A GLOBAL VARIABLE.')
  fprintf('Type ''dbcont'' to run startup.m and continue\n or ''dbquit'' to fix this yourself ...\n');
  keyboard;
  startup
end

% LOAD THE PARAM SPREADSHEET AND SAVE LOCAL VARIABLES FROM THE SETTINGS
params = read_param_xls(ct_filename_param(settings.paramFn));
settings.radarName = params(1).radar_name;
settings.seasonName = params(1).season_name;

% Join sys_segments and sys_point_paths to calculate the number of point paths in each segment
try
  query1 = sprintf('SELECT seg.id, seg.name, COUNT(pp.id) FROM %s_segments seg JOIN %s_point_paths pp ON seg.id=pp.segment_id WHERE seg.season_id = (SELECT id FROM %s_seasons WHERE name = ''%s'') GROUP BY seg.id;', settings.sysName, settings.sysName, settings.sysName, settings.seasonName);
  [~,data1] = opsQuery(query1);
  segment_ID = [data1{1,:}]; % If each cell contains the same type of data, you can create a single variable by applying the array concatenation operator, [], to the comma-separated list.  % Two single quote instead of one double quote
  seg = sprintf('%d,',segment_ID(:));
  seg = seg(1:end-1);
catch ME
  ME.getReport();
end

% Join sys_layer_points and sys_point_paths to calculate the number of picked points
try
  query2 = sprintf('SELECT pp.segment_id, lp.layer_id, COUNT (lp.id) FROM %s_layer_points lp JOIN %s_point_paths pp ON lp.point_path_id=pp.id WHERE pp.segment_id IN (%s) AND lp.layer_id IN (1,2) GROUP BY pp.segment_id, lp.layer_id;', settings.sysName, settings.sysName, seg);
  [~,data2] = opsQuery(query2);
catch ME
  ME.getReport();
end

% Report the status of segment insertion and picking. The report includes
% four fields: seg_Name (segment name), inserted (flag whether this segment
% has been inserted or not), surfacePickedPercent (the percentage of
% surface layer has been picked), bottomPickedPercent (the percentage of
% bottom layer has been picked)
try
  report = [];
  number = 1;
  for param_idx = 1:length(params)
    
    param = params(param_idx);
    if param.cmd.generic == 1
      report(number).seg_Name = param.day_seg;
      
      for i = 1:size(data1,2)
        
        TF = strcmp(param.day_seg, data1{2,i});
        if TF
          report(number).inserted = 'Yes';
          for j = 1:size(data2,2)
            if (data2{1,j} == data1{1,i})&& (data2{2,j} == 1)
              report(number).surfacePickedPercent = double(data2{3,j})/double(data1{3,i})*100.0;
              report(number).bottomPickedPercent = double(data2{3,j+1})/double(data1{3,i})*100.0;
              break;
            end
          end
          break;
        end
        
        if i == size(data1,2)
          report(number).inserted = 'No';
          report(number).surfacePickedPercent = 0;
          report(number).bottomPickedPercent = 0;
        end
        
      end
      
      number = number + 1;
    end
    
  end
catch ME
  fprintf('\n');
  warning(sprintf('%s at line %d in file %s.',ME.message,ME.stack(1).line,ME.stack(1).name));
end

% Remove all variables from the workspace except for the variable report
clearvars -except report;