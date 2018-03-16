% script runOpsShiftGpsTime
%
% Function for running opsShiftGpsTime
%
% Shift the GPS time of layer points in the database.

params = read_param_xls(ct_filename_param('rds_param_2010_Greenland_DC8.xls'),'20100414_02','post');
% Layer to shift
lyr_name = 'surface';
% Offset to add to the GPS time (sec)
offset = 0;
% Quality to override (NaN: do not change, 1: good, 2: medium, 3: poor)
quality = 2;

params(1).cmd.generic = 1;

%% Automated Section
% =====================================================================
% Create param structure array
% =====================================================================

%% Process each of the segments
% =====================================================================
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  sys = ct_output_dir(param.radar_name);
  
  shift_param = [];
  shift_param.properties.location = param.post.ops.location;
  shift_param.properties.season = param.season_name;
  shift_param.properties.lyr_name = lyr_name;
  shift_param.properties.segment = param.day_seg;
  shift_param.properties.offset = offset;
  shift_param.properties.quality = quality;
  fprintf('GPS time shift %s by %.1f sec (%s)\n', param.day_seg, offset, datestr(now));
  opsShiftGpsTime(sys,shift_param);
  
end
fprintf('Completed (%s)\n', datestr(now));
