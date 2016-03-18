% script runOpsShiftTwtt
%
% Function for running opsShiftTwtt
%
% Shift the TWTT of layer points in the database.

params = read_param_xls(ct_filename_param('rds_param_2010_Greenland_DC8.xls'),[],'post');
% Layer to offset
lyr_name = 'bottom';
% Offset to add to the twtt (sec)
offset = NaN;

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
  
  if ~isnan(offset)
    shifttwtt_param = [];
    shifttwtt_param.properties.location = param.post.ops.location;
    shifttwtt_param.properties.season = param.season_name;
    shifttwtt_param.properties.lyr_name = lyr_name;
    shifttwtt_param.properties.segment = param.day_seg;
    shifttwtt_param.properties.offset = offset;
    fprintf('Shifting %s by %.1f ns\n', param.day_seg, offset*1e9);
    opsShiftTwtt(sys,shifttwtt_param);
    
  else
    delete_param = [];
    delete_param.properties.location = param.post.ops.location;
    delete_param.properties.season = param.season_name;
    delete_param.properties.lyr_name = lyr_name;
    delete_param.properties.segment = param.day_seg;
    % This is a hack to select "all" points
    %   Decimal offsets are required so that python treats as double rather
    %   than integer.
    delete_param.properties.start_gps_time = 0.1;
    delete_param.properties.stop_gps_time = 1e14+0.1;
    delete_param.properties.min_twtt = -10.1;
    delete_param.properties.max_twtt = 10.1;
    delete_param.properties.lyr_name = lyr_name;
    warning('Deleting %s. Type dbcont to continue.\n', param.day_seg);
    keyboard
    opsDeleteLayerPoints(sys,delete_param);
  end
  
end
