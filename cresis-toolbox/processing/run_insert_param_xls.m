% script run_insert_param_xls
%
% Example to run insert_param_xls.m
%
% Author: Jordan Sprick

param_fn = ct_filename_param('rds_param_2011_Greenland_P3.xls');

% Make a test copy of the param spreadsheet
[param_fn_dir,param_fn_name,param_fn_ext] = fileparts(param_fn);
new_param_fn = fullfile(param_fn_dir,[param_fn_name '_test' param_fn_ext]);
copyfile(param_fn,new_param_fn);

% Read in the original parameter spreadsheet
param_original = read_param_xls(param_fn);

% Modify the structure
param_new = param_original;

% Add a new radar waveform to each of the original segments
for p = 1:length(param_original)
  param_new(p).radar.wfs(end+1).Tpd = 100e-6;
end

% Modify the param spreadsheet
insert_param_xls(new_param_fn,param_new,'radar');

% Read in the modified param spreadsheet
param_read_new = read_param_xls(new_param_fn);

% Add a new segment to the last day
last_day_seg = param_new(end).day_seg;
new_day_seg = sprintf('%s%02d',last_day_seg(1:end-2),str2num(last_day_seg(end-1:end))+1);
param_new(end+1).day_seg = new_day_seg;

% Modify the param spreadsheet
insert_param_xls(new_param_fn,param_new);

% Read in the modified param spreadsheet
param_read_new2 = read_param_xls(new_param_fn);

return;
