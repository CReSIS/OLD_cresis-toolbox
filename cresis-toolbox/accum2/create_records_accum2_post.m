% script create_records_accum2_post
%
% This is a script which is called by create_records_accum2
% Author: John Paden, Anthony Hoch, Jilu Li, Logan Smith

dbstack_info = dbstack;
if ~exist('param','var') || isempty(param) || length(dbstack_info) == 1
  % =====================================================================
  % Debug Setup
  % =====================================================================
  new_param = read_param_xls(ct_filename_param('accum_param_2013_Antarctica_Ground.xls'),'20140110_01');

  fn = ct_filename_tmp(new_param,new_param.records.records_fn,'records','workspace');
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

fprintf('Running %s gps sync (%s)\n', param.day_seg, datestr(now));

% ===================================================================
% ===================================================================
% Correlate GPS with radar data
% ===================================================================
% ===================================================================

% Apply a time offset (useful when radar system does not have GPS time
% lock).  For ideal radar operation time_offset is zero.
if param.vectors.gps.time_offset ~= 0
  warning('ACCUM2 usually has 0 time offset, but is set to %g', param.vectors.gps.time_offset);
end
records.radar_time = records.radar_time + param.vectors.gps.time_offset;

%% Look for Bad Records
% State machine has 3 states:
% 0: init --> 1 when good record is found
% 1: last record was good --> 2 when bad jump is found
% 2: last record had bad jump --> 1 when good jump is found
%   Side effect: when transitioning back to 1, bad records are flagged
fstate = 0;
bad_mask = zeros(size(records.radar_time));
first_good_rec = NaN;
num_coh_ave_start = median(records.num_coh_ave(1:7));
num_coh_ave_end = median(records.num_coh_ave(end-6:end));
epri = length(records.wfs{1}) * medfilt1(records.num_coh_ave,7)/param.radar.prf;
records.num_coh_ave(1:4) = num_coh_ave_start;
records.num_coh_ave(end-3:end) = num_coh_ave_end;
tepri = diff(records.radar_time);
good_jump = abs(tepri - epri(1:end-1)) < epri(1:end-1)/10;
for cur_rec = 1:length(records.radar_time)-1
  if good_jump(cur_rec)
    if fstate ~= 1
      if fstate == 0
        % Found first good record, flag all records before this as bad
        fstate = 1;
        first_good_rec = cur_rec;
        if length(1:cur_rec-1) > 0
          warning('Found bad jumps (%d dropped records)\n', cur_rec-1);
        end
        bad_mask(1:cur_rec-1) = 1;
      elseif fstate == 2
        % Found next good record, flag bad records
        if cur_rec - last_good_rec - 1 > 0
          warning('Found bad jumps (%d dropped records)\n', cur_rec - last_good_rec - 1);
        end
        bad_mask(last_good_rec+1:cur_rec-1) = 1;
        fstate = 1;
      end
    end
  else
    if fstate == 1
      % Search for next good record
      fstate = 2;
      last_good_rec = cur_rec;
    end
  end
end
if fstate == 2
  bad_mask(last_good_rec+1:end) = 1;
end
if fstate == 0
  warning('No good records found, debug reason\n');
  % No good records
  keyboard
  error('No good records found');
end
if any(bad_mask)
  warning('Found %d bad records', sum(bad_mask));
end

%% Drop Bad Records
records.comp_time = records.comp_time(~bad_mask);
records.radar_time = records.radar_time(~bad_mask);
records.radar_time_1pps = records.radar_time_1pps(~bad_mask);
records.range_gate_start = records.range_gate_start(~bad_mask);
records.range_gate_duration = records.range_gate_duration(~bad_mask);
records.num_coh_ave = records.num_coh_ave(~bad_mask);
for adc_idx = 1:length(records.file_idx)
  records.file_idx{adc_idx} = records.file_idx{adc_idx}(~bad_mask);
  records.offset{adc_idx} = records.offset{adc_idx}(~bad_mask);
end

%% Correct First Record In File Offsets (for dropping of bad records)
for file_idx = 1:length(records.file_rec_offset)
  % The records.file_rec_offset vector contains the record position of the
  % first record in each file.  We need to adjust this for the dropping of
  % bad records
  num_bad_before_this_file = sum(bad_mask(1:records.file_rec_offset(file_idx)-1));
  records.file_rec_offset(file_idx) ...
    = records.file_rec_offset(file_idx) - num_bad_before_this_file;
  if records.file_rec_offset(file_idx) > length(records.radar_time)
    warning('File %d is no longer needed (all bad records)', file_idx);
  end
end

if param.records.gps.en
  records = sync_radar_to_gps(param,records,records.radar_time,records.comp_time);
else
  % GPS file does not exist, so fill in with NaN
  warning('GPS file %s does not exist (writing radar UTC time as GPS time!)', gps_fn);
  records.lat = NaN*zeros(size(records.comp_time));
  records.lon = NaN*zeros(size(records.comp_time));
  records.elev = NaN*zeros(size(records.comp_time));
  records.gps_time = NaN*zeros(size(records.comp_time));
  records.roll = NaN*zeros(size(records.comp_time));
  records.pitch = NaN*zeros(size(records.comp_time));
  records.heading = NaN*zeros(size(records.comp_time));
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

board_idx = 1;
records.relative_filename = [];
records.relative_rec_num = [];
tmp = records.offset{board_idx};
records.offset = [];
unique_file_idxs = unique(records.file_idx{board_idx});
for file_idx = 1:length(unique_file_idxs)
  records.relative_rec_num{board_idx}(file_idx) = find(records.file_idx{board_idx} == unique_file_idxs(file_idx),1);
  [fn_dir fn_name fn_ext] = fileparts(records.filenames{board_idx}{unique_file_idxs(file_idx)});
  records.relative_filename{board_idx}{file_idx} = [fn_name fn_ext];
end
records = rmfield(records,'file_idx');
records = rmfield(records,'filenames');
records = rmfield(records,'file_rec_offset');
records.offset(board_idx,:) = tmp;

records.radar_name = param.radar_name;
records.ver = 3;
%  records.raw.epri(1 .. Nx)
%  records.raw.seconds(1 .. Nx)
%  records.raw.fraction(1 .. Nx)
records.raw.comp_time = records.comp_time;
records.raw.radar_time = records.radar_time;
records.raw.radar_time_1pps = records.radar_time_1pps;
records = rmfield(records,'comp_time');
records = rmfield(records,'radar_time');
records = rmfield(records,'radar_time_1pps');

records.settings = [];

records.settings.wfs_records = 1;
records.settings.wfs = records.wfs{1};
% for wf = 1:8 % HACK
%   records.settings.wfs(wf) = records.settings.wfs(1);
% end
for wf = 1:length(records.settings.wfs)
  records.settings.wfs(wf).bit_shifts = log(records.settings.wfs(wf).presums)/log(2);
end
records.settings.range_gate_start = records.range_gate_start;
records.settings.range_gate_duration = records.range_gate_duration;
records.settings.trigger_delay = records.trigger_delay(1) * ones(size(records.range_gate_start));
records.settings.num_coh_ave = records.num_coh_ave;
records = rmfield(records,'wfs_file');
records = rmfield(records,'wfs');
records = rmfield(records,'range_gate_start');
records = rmfield(records,'range_gate_duration');
records = rmfield(records,'trigger_delay');
records = rmfield(records,'num_coh_ave');

records.notes = '';
records.param_records = param;

fprintf('Saving records file %s (%s)\n',records_fn,datestr(now));
save(records_fn,'-v7.3','-struct','records'); % Handle large file sizes, so use v7.3

% =====================================================================
% Create record aux files for faster loading times
% =====================================================================

fprintf('Creating auxiliary records files %s (%s)\n',records_fn,datestr(now));
create_records_aux_files(records_fn);

fprintf('Done (%s)\n\n', datestr(now));
