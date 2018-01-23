% script create_records_acords_post
%
% This is a script which is called by create_records_mcrds 
% Author: John Paden

% ===================================================================
% ===================================================================
% Correlate GPS with radar data
% ===================================================================
% ===================================================================
num_adcs = length(param.records.file.adcs);
for adc_idx = 1:num_adcs
  % ===================================================================
  % Arrange the output variables and synchronize with the GPS
  if adc_idx == 1
    records.epri = uint32(hdrs.epri{adc_idx});
    records.time = hdrs.time{adc_idx};
    records.time        = datenum(datevec(records.time/(24*60*60) + datenum(1970,1,1,0,0,0)));     
    records.time        = records.time + gps.time_offset/(24*60*60);                                        
    [year month day hour minute sec] = datevec(records.time);
    records.time        = hour*3600 + minute*60 + sec;                  % UTC seconds of day   
    records.time = records.time';
    records.ver = hdrs.ver{adc_idx}(1);
    records.wfs = wfs;
    
    % Check for seconds of day roll over and unwrap
    wrap_idxs = find(diff(records.time) < -86300);
    for wrap_idx = wrap_idxs
      records.time(wrap_idx+1:end) = records.time(wrap_idx+1:end) + 86400;
    end
    
    if isfield(param.records.gps,'time_offset') && ~isempty(param.records.gps.time_offset)
      % Apply a time offset (useful when radar system does not have GPS time
      % lock).  For ideal radar operation time_offset is zero.
      records.time = records.time + param.records.gps.time_offset;
    end
    
    if param.records.gps.en
      % Synchronize times to get positions and absolute time
      records.lat = double(interp1(UTC_sod,gps.lat,records.time));
      records.lon = double(mod(interp1(UTC_sod,unwrap(gps.lon/180*pi),records.time)*180/pi+180, 360)-180);
      records.elev = double(interp1(UTC_sod,gps.elev,records.time));
      records.roll = double(interp1(UTC_sod,gps.roll,records.time));
      records.pitch = double(interp1(UTC_sod,gps.pitch,records.time));
      records.heading = double(mod(interp1(UTC_sod,unwrap(gps.heading),records.time)+pi,2*pi)-pi);
      records.gps_time = interp1(UTC_sod,gps.gps_time,records.time);
    end
  end
  records.file_idx{adc_idx} = hdrs.file_idx{adc_idx};
  records.offset{adc_idx} = hdrs.offset{adc_idx};
end

records.gps_source = gps.gps_source;

% ===================================================================
% ===================================================================
% Calculate final variables and save output
% ===================================================================
% ===================================================================
fprintf('Create variables (%.1f sec)\n', toc);

records.filenames = filenames;

records.eprf_computed = 1/mean(diff(hdrs.time{1}));
records.prf_computed = records.eprf_computed * num_presum;

% Convert absolute filenames into relative filenames
for adc_idx = 1:length(records.filenames)
  for file_idx = 1:length(records.filenames{adc_idx})
    [path name ext] = fileparts(records.filenames{adc_idx}{file_idx});
    records.filenames{adc_idx}{file_idx} = [name ext];
  end
end

% Make sure output directory exists
[out_dir out_name] = fileparts(ct_filename_support(param,'','records'));
if ~exist(out_dir,'dir')
  mkdir(out_dir);
end

fprintf('Save variables (%.1f sec)\n', toc);
param_records = param.records;
param_radar = param.radar;
records_fn = ct_filename_support(param,'','records');
save(records_fn,'records','param_records','param_radar');
create_records_aux_files(records_fn);

% Make sure output directory exists
[out_dir out_name] = fileparts(ct_filename_support(param,param.frames.frames_fn,'frames'));
if ~exist(out_dir,'dir')
  mkdir(out_dir);
end

return;
