function check_records(param)
% check_records(param)
%
% Checks the fields in the records files for potential errors.
%
% param = parameter structure indicating which radar and season to load.
%   This causes the standard parameter spreadsheet filename to be created
%   and loaded.  Therefore it will not work if you don't have the parameter
%   spreadsheet in the standard location with the standard name.
%   ALTERNATE MODE: pass in a string containing the filename of a specific
%   record to load
%
% Example:
% param = [];
% param.radar_name = 'mcrds';
% param.season_name = '2009_Antarctica_TO';
% check_records(param)
%
% param = [];
% param.radar_name = 'mcords2';
% param.season_name = '2011_Greenland_P3';
% check_records(param)
%
% param = [];
% param.radar_name = 'snow3';
% param.season_name = '2013_Antarctica_Basler';
% check_records(param)

if ischar(param)
  % param is a string containing the filename
  records = load(param,'param_records');
  if ~isfield(records,'param_records')
    error('This mode only supported for new records files with param_records field.');
  end
  check_records_support_func(param,records.param_records,-inf);
  
elseif isstruct(param)
  % param is a struct indicating which radar/season to check (checks all segments)
  [output_dir,radar_type] = ct_output_dir(param.radar_name);
  
  param_fn = ct_filename_param(sprintf('%s_param_%s.xls', output_dir, param.season_name));
  params = read_param_xls(param_fn);
  
  old_time = -inf;
  for param_idx = 1:length(params)
    param = params(param_idx);
    
    records_fn = ct_filename_support(param,'','records');
    
    if ~isempty(regexpi(param.cmd.notes,'do not process'))
      continue;
    end
    
    % DEBUG OPTION TO JUST CHECK RECORDS BASED ON GENERIC COLUMN IN SPREADSHEET
    % Uses the generic column of the parameter spreadsheet to determine which segments
    % to check.
    %if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    %  continue;
    %end

    fprintf('%d of %d: Checking records %s\n', param_idx, length(params), records_fn);
    
    old_time = check_records_support_func(records_fn,param,old_time);
  end
  
else
  error('Invalid argument')
end

end

function old_time = check_records_support_func(records_fn,param,old_time)
% check_records_support_func(records_fn,param,old_time)
%
% Support function which does the actual checking of the records

if ~exist(records_fn,'file')
  warning('Missing records file %s\n', records_fn);
  return;
end

records = load(records_fn);

if isfield(records,'records')
  warning('Old records file format, skipping');
  return;
end

if isempty(regexpi(records.gps_source,'final'))
  warning('GPS source is %s and not final\n', records.gps_source);
end

if strcmpi(records.radar_name,'mcrds')
  if length(records.raw.wfs) > 1
    warning('Waveform headers change during records');
  end
  
  if isfield(records,'adc_phase_corr_deg')
    fprintf('  This record has adc_phase_corr_deg field\n');
  end
  
  radar_clock_delta = records.raw.radar_time - records.raw.comp_time;
  radar_clock_delta = detrend(radar_clock_delta);
  radar_clock_delta = radar_clock_delta - radar_clock_delta(1);
  
  if any(abs(radar_clock_delta) > 0.5)
    warning('Radar and computer clocks show too large delta');
    plot(radar_clock_delta)
    keyboard;
  end
  
  settings_inconsistent = false;
  for settings_idx = 1:length(records.raw.wfs)
    for wf = 1:length(records.raw.wfs(settings_idx).Waveform)
      if compare_structs(records.raw.wfs(settings_idx).Waveform(wf), records.raw.wfs(1).Waveform(wf),1)
        settings_inconsistent = true;
      end
      if records.raw.wfs(settings_idx).Waveform(wf).SampleDelay(1) ~= records.raw.wfs(1).Waveform(wf).SampleDelay(1)
        warning('SampleDelay not consistent');
      end
      if records.raw.wfs(settings_idx).Waveform(wf).NumberSamples(1) ~= records.raw.wfs(1).Waveform(wf).NumberSamples(1)
        warning('NumberSamples not consistent');
      end
      if param.radar.wfs(wf).Tpd ~= records.raw.wfs(settings_idx).Waveform(wf).PulseDuration
        warning('%d: Pulse duration does not match',settings_idx);
      end
      if param.radar.wfs(wf).f0 ~= records.raw.wfs(settings_idx).Waveform(wf).StartFrequency
        warning('%d: Start frequency does not match',settings_idx);
      end
      if param.radar.wfs(wf).f1 ~= records.raw.wfs(settings_idx).Waveform(wf).StopFrequency
        warning('%d: Stop frequency does not match',settings_idx);
      end
      if 20*log10(param.radar.wfs(wf).adc_gains(1)) ~= 70.0-records.raw.wfs(settings_idx).Waveform(wf).RxAttenuation
        warning('%d: Receiver attenuation does not match',settings_idx);
      end
    end
  end
  if settings_inconsistent
    records.raw.wfs_file
  end
  
end

along_track = geodetic_to_along_track(records.lat,records.lon,records.elev);

vel = diff(along_track) ./ diff(records.gps_time);
if any(vel < 25 | vel > 300)
  if all(vel < 6)
    warning('Small or large velocity. All below 6 m/s.');
  else
    warning('Small or large velocity. Some above 6 m/s.');
  end
end

if records.lat >= 90 | records.lat <= -90
  warning('lat out of bounds');
end

if records.lon >= 360 | records.lon <= -360
  warning('lon out of bounds');
end

if records.elev >= 40000 | records.elev <= -10000
  warning('elev out of bounds');
end

if records.roll >= 100/180*pi | records.roll <= -100/180*pi
  warning('roll out of bounds');
end

if records.pitch >= pi/2 | records.pitch <= -pi/2
  warning('pitch out of bounds');
end

if records.heading >= 2*pi | records.heading <= -2*pi
  warning('heading out of bounds');
end

if old_time > records.gps_time(1)
  warning('records out of order');
end

if any(diff(records.gps_time) < 0)
  warning('time not monotonically increasing')
end

if any(isnan(records.lat) | isnan(records.lon) | isnan(records.elev)  | isnan(records.roll)  | isnan(records.pitch)  | isnan(records.heading))
  warning('NaN in record');
end

diff_time = diff(records.gps_time);
if any(diff_time > 1000 / median(vel))
  jump_idxs = find(diff_time > 1000 / median(vel));
  warning('Time gap greater than 1 km in record (assuming %f m/s)', median(vel));
  fprintf('Time gaps in seconds: '); fprintf('%f\t', diff_time(jump_idxs)); fprintf('\n');
  for jump_idx = jump_idxs
    fprintf('  Gap is at file index %d to %d\n', find(records.relative_rec_num{1} > jump_idx,1), ...
      find(records.relative_rec_num{1} > jump_idx+1,1));
  end
end

old_time = records.gps_time(end);

end
