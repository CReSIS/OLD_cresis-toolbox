% script create_records_mcrds_sync
%
% This is a script which is called by create_records_mcrds
% Author: John Paden, Anthony Hoch, Jilu Li

dbstack_info = dbstack;
if ~exist('param','var') || isempty(param) || length(dbstack_info) == 1
  % =====================================================================
  % Debug Setup
  % =====================================================================
  new_param = read_param_xls(ct_filename_param('rds_param_2009_Greenland_TO.xls'),'20090402_01');

  fn = ct_filename_ct_tmp(new_param,'','records','workspace');
  fn = [fn '.mat'];
  fprintf('Loading workspace %s (%s)\n', fn, datestr(now));
  if exist(fn,'file')
    load(fn);
  else
    error('Temporary records file does not exist');
  end
  param = new_param;
  clear new_param;
  
  clear('param_override');
  param_override.sched.type = 'no scheduler';
  
  % Input checking
  if ~exist('param','var')
    error('A struct array of parameters must be passed in\n');
  end
  global gRadar;
  if exist('param_override','var')
    param_override = merge_structs(gRadar,param_override);
  else
    param_override = gRadar;
  end
  param = merge_structs(param, param_override);
end

fprintf('Running %s correction and gps sync (%s)\n', param.day_seg, datestr(now));

% ===================================================================
% ===================================================================
% Correlate GPS with radar data
% ===================================================================
% ===================================================================

% Apply a time offset (useful when radar system does not have GPS time
% lock).  For ideal radar operation time_offset is zero.
if param.vectors.gps.time_offset ~= 0
  warning('MCRDS usually has 0 time offset, but is set to %g', param.vectors.gps.time_offset);
end
hdrs.radar_time = hdrs.radar_time + param.vectors.gps.time_offset;

records = [];
if param.records.gps.en
  records = sync_radar_to_gps(param,records,hdrs.radar_time,hdrs.comp_time); 
else
  % GPS file does not exist, so fill in with NaN
  warning('GPS file %s does not exist (writing radar UTC time as GPS time!)', gps_fn);
  records.lat = NaN*zeros(size(records.radar_time));
  records.lon = NaN*zeros(size(records.radar_time));
  records.elev = NaN*zeros(size(records.radar_time));
  records.gps_time = NaN*zeros(size(records.radar_time));
  records.roll = NaN*zeros(size(records.radar_time));
  records.pitch = NaN*zeros(size(records.radar_time));
  records.heading = NaN*zeros(size(records.radar_time));
  records.gps_source = 'NA';
end

% =====================================================================
% Save records files
% =====================================================================

% Save concatenated files in records directories after time fix
records_fn = ct_filename_support(param,'','records');
[records_fn_dir records_fn_name] = fileparts(records_fn);
if ~exist(records_fn_dir,'dir')
  fprintf('Output directory %s does not exist, creating...\n', records_fn_dir);
  mkdir(records_fn_dir);
end

% Standard Fields
% records
%  .lat
%  .lon
%  .elev
%  .roll
%  .pitch
%  .heading
%  .gps_time
%  .gps_source
records.surface = NaN*zeros(size(records.lat));
records.relative_filename{1} = hdrs.filenames;
records.relative_rec_num{1} = hdrs.file_rec_offset;

% Create records.offset variable
rec_data_size = hdrs.wfs(1).IndexRecordStop - hdrs.wfs(1).IndexRecordStart + 1;
rec_data_size = sum(rec_data_size(:));
header_size = 484 + 176*hdrs.wfs(1).NumberWaveforms;
rec_header_size = 24;
rec_size = (rec_header_size+2*rec_data_size);
records.offset = [];
for file_idx = 1:length(hdrs.file_rec_offset)-1
  num_records = hdrs.file_rec_offset(file_idx+1) - hdrs.file_rec_offset(file_idx);
  records.offset = [records.offset (header_size+rec_header_size+rec_size*(0 : num_records-1))];
end
records.offset = [records.offset (header_size+rec_header_size+rec_size*(0 : length(records.lat)-length(records.offset)-1))];

records.radar_name = param.radar_name;
records.ver = 3;
records.raw.comp_time = hdrs.comp_time;
records.raw.radar_time = hdrs.radar_time;
records.raw.wfs = hdrs.wfs;
records.raw.wfs_file = hdrs.wfs_file;

records.settings = [];

% Copy waveform structs over into standard format
records.settings.wfs_records = records.relative_rec_num{1}(hdrs.wfs_file);
for wfs_idx = 1:length(hdrs.wfs)
  for wf = 1:length(hdrs.wfs(wfs_idx).Waveform)
    records.settings.wfs(wfs_idx).wfs(wf).presums = hdrs.wfs(wfs_idx).PreSum;
    records.settings.wfs(wfs_idx).wfs(wf).bit_shifts ...
      = floor(log((records.settings.wfs(wfs_idx).wfs(wf).presums-1)/16)/log(2)) + 1;
    if records.settings.wfs(wfs_idx).wfs(wf).bit_shifts < 0
      records.settings.wfs(wfs_idx).wfs(wf).bit_shifts = 0;
    end
    records.settings.wfs(wfs_idx).wfs(wf).fs = hdrs.wfs(wfs_idx).SampleFrequency;
    records.settings.wfs(wfs_idx).wfs(wf).f0 = hdrs.wfs(wfs_idx).Waveform(wf).StartFrequency;
    records.settings.wfs(wfs_idx).wfs(wf).f1 = hdrs.wfs(wfs_idx).Waveform(wf).StopFrequency;
    records.settings.wfs(wfs_idx).wfs(wf).Tpd = hdrs.wfs(wfs_idx).Waveform(wf).PulseDuration;
    % Currently assumes that all ADCs have the same sample delay
    if any(hdrs.wfs(wfs_idx).Waveform(wf).SampleDelay ~= hdrs.wfs(wfs_idx).Waveform(wf).SampleDelay(1))
      error('Code does not handle different SampleDelay on each ADC');
    end
    records.settings.wfs(wfs_idx).wfs(wf).Tadc = hdrs.wfs(wfs_idx).Waveform(wf).SampleDelay(1);
    records.settings.wfs(wfs_idx).wfs(wf).start_idx = hdrs.wfs(wfs_idx).Waveform(wf).SampleDelay(1)*hdrs.wfs(wfs_idx).SampleFrequency;
    records.settings.wfs(wfs_idx).wfs(wf).num_sam = hdrs.wfs(wfs_idx).Waveform(wf).NumberSamples(1);
    records.settings.wfs(wfs_idx).wfs(wf).adc_gains = 10.^((72 - hdrs.wfs(wfs_idx).Waveform(wf).RxAttenuation)/20);
    records.settings.wfs(wfs_idx).wfs(wf).adc_en = hdrs.wfs(wfs_idx).Waveform(wf).RecordEnable;
  end
end

records.notes = '';
records.param_records = param;

fprintf('Saving records file %s (%s)\n',records_fn,datestr(now));
save(records_fn,'-v6','-struct','records');

% =====================================================================
% Create record aux files for faster loading times
% =====================================================================

fprintf('Creating auxiliary records files %s (%s)\n',records_fn,datestr(now));
create_records_aux_files(records_fn);

fprintf('Done (%s)\n\n', datestr(now));

return
