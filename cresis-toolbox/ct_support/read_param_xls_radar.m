function [params] = read_param_xls_radar(param_fn, day_seg_filter)
% [params] = read_param_xls_radar(param_fn, day_seg_filter)
%
% Support function for read_param_xls (for acords, mcrds, mcords,
% mcords2 radars).
%
% Author: Brady Maasen, John Paden
%
% See also: ct_set_params, master, read_param_xls
%
% See also for spreadsheet cell loading:
%  read_param_xls_boolean.m, read_param_xls_general.m,
%  read_param_xls_text.m
%  
% See also for worksheet loading:
%  read_param_xls_generic.m, read_param_xls_radar.m: 
%
% See also for printing out spreadsheet to stdout:
%  read_param_xls_print, read_param_xls_print_headers.m

cell_boolean = @read_param_xls_boolean;
cell_text = @read_param_xls_text;
cell_read = @read_param_xls_general;

warning('off','MATLAB:xlsread:Mode');

%% Create Command Parameters
% =======================================================================
fprintf('Reading command/cmd of xls file: %s\n', param_fn);

warning off MATLAB:xlsfinfo:ActiveX
[status, sheets] = xlsfinfo(param_fn);
warning on MATLAB:xlsfinfo:ActiveX
cmd_sheet_name = 'cmd';
sheet_idx = strmatch(cmd_sheet_name,sheets,'exact');
if isempty(sheet_idx)
  cmd_sheet_name = 'command';
  sheet_idx = strmatch(cmd_sheet_name,sheets,'exact');
  if isempty(sheet_idx)
    error('  Command/cmd sheet not found.');
  end
end

[num, txt] = xlsread(param_fn,cmd_sheet_name,'','basic');

param_file_version        = cell_text(1,2,num,txt);
radar_name                = cell_text(2,2,num,txt);
season_name               = cell_text(3,2,num,txt);
sw_version                = current_software_version;
if ispc
  user_name = getenv('USERNAME');
else
  [~,user_name] = system('whoami');
  user_name = user_name(1:end-1);
end

% Interpret version field: text in row 1, column 2
version = str2double(param_file_version);
if isnan(version)
  error('Could not find version field text in row 1, column 2 of the first worksheet.');
end
if version < 4
  num_header_rows = 5;
  rows = max(size(num,1), size(txt,1)) - num_header_rows;
  
  for idx = 1:rows
    params(idx).fn                                = param_fn;
    params(idx).user_name                         = user_name;
    params(idx).radar_name                        = radar_name;
    params(idx).season_name                       = season_name;
    row = idx + num_header_rows;
    params(idx).day_seg                           = sprintf('%08.0f_%02.0f',num(row,1),num(row,2));
    col = 3;
    params(idx).cmd.frms                          = cell_read(row,col,num,txt); col = col + 1;
    params(idx).cmd.vectors                       = cell_boolean(row,col,num,txt); col = col + 1;
    params(idx).cmd.records                       = cell_boolean(row,col,num,txt); col = col + 1;
    params(idx).cmd.frames                        = cell_boolean(row,col,num,txt); col = col + 1;
    params(idx).cmd.qlook                         = cell_boolean(row,col,num,txt); col = col + 1;
    params(idx).cmd.sar                           = cell_boolean(row,col,num,txt); col = col + 1;
    params(idx).cmd.array                         = cell_boolean(row,col,num,txt); col = col + 1;
    params(idx).cmd.generic                       = cell_read(row,col,num,txt); col = col + 1;
    if isempty(params(idx).cmd.generic)
      params(idx).cmd.generic = 0;
    end
    params(idx).cmd.mission_names                 = cell_text(row,col,num,txt); col = col + 1;
    params(idx).cmd.notes                         = cell_text(row,col,num,txt); col = col + 1;
    params(idx).sw_version                        = sw_version;
    params(idx).param_file_version                = param_file_version;
  end

else
  sheet_name = 'cmd';
  [params] = read_param_xls_generic(param_fn,sheet_name,[]);
  for idx = 1:length(params)
    if isempty(params(idx).cmd.generic)
      params(idx).cmd.generic = 0;
    end
    params(idx).fn                                = param_fn;
    params(idx).user_name                         = user_name;
    params(idx).radar_name                        = radar_name;
    params(idx).season_name                       = season_name;
    params(idx).sw_version                        = sw_version;
    params(idx).param_file_version                = param_file_version;
  end
end

%% Create Records, Frames Parameters
% =======================================================================
sheet_name = 'records';

[params] = read_param_xls_generic(param_fn,sheet_name,params);

%% qlook parameters
% =======================================================================
sheet_name = 'qlook';

[params] = read_param_xls_generic(param_fn,sheet_name,params);

%% SAR parameters
% =======================================================================
sheet_name = 'sar';

[params] = read_param_xls_generic(param_fn,sheet_name,params);

%% Combine waveforms and channel parameters
% =======================================================================
sheet_name = 'array';

[params] = read_param_xls_generic(param_fn,sheet_name,params);

%% Radar configuration parameters
% =======================================================================
sheet_name = 'radar';

[params] = read_param_xls_generic(param_fn,sheet_name,params);

%% Data posting parameters
% =======================================================================
sheet_name = 'post';

[params] = read_param_xls_generic(param_fn,sheet_name,params);
