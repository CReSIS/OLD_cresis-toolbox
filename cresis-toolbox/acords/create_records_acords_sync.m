% create_records_fmcw_acords_sync
%
% Script called from create_records_fmcw_acords
%
% To run this in debug mode, you need to set the debug setup section in
% this file and then just run this as a script (must be run from ">>"
% prompt and not a debug "K>>" prompt).
%
% Author: John Paden

dbstack_info = dbstack;
if ~exist('param','var') || isempty(param) || length(dbstack_info) == 1
  % =====================================================================
  % Debug Setup
  % =====================================================================
  
  new_param = read_param_xls(ct_filename_param('rds_param_2005_Greenland_TO.xls'),'20050505_01');

  fn = ct_filename_tmp(new_param,'','records','workspace');
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

% ======================================================================

if ~isfield(param.records,'manual_time_correct') || isempty(param.records.manual_time_correct)
  param.records.manual_time_correct = 0;
end

radar_time_notes = '';

%% Create hdr.comp_time vector
if param.records.file_version == 101
  hdr.comp_time = double(hdr.seconds) + 2*double(hdr.fraction)/param.radar.fs;
else
  hdr.comp_time = double(hdr.seconds) + double(hdr.fraction)/param.radar.fs;
end

% ===================================================================
%% Correlate GPS with radar data
% ===================================================================
fprintf('Loading GPS data (%s)\n', datestr(now));

if param.records.gps.en
  records = sync_radar_to_gps(param,records,hdr.comp_time);
  
else
  records.lat = NaN*zeros(size(hdr.comp_time));
  records.lon = NaN*zeros(size(hdr.comp_time));
  records.elev = NaN*zeros(size(hdr.comp_time));
  records.gps_time = NaN*zeros(size(hdr.comp_time));
  records.roll = NaN*zeros(size(hdr.comp_time));
  records.pitch = NaN*zeros(size(hdr.comp_time));
  records.heading = NaN*zeros(size(hdr.comp_time));
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
%  .relative_filename{board_idx == 1}{1 .. Nf(1)}
%  .relative_rec_num{board_idx == 1}(1 .. Nf(1))
records.offset = hdr.offset;
records.radar_name = param.radar_name;
records.ver = 3;
% if param.records.file_version ~= 101
%   records.raw.epri = orig_hdr.epri;
% end
% records.raw.seconds = orig_hdr.seconds;
% records.raw.fraction = orig_hdr.fraction;
records.raw.seconds = hdr.seconds;
records.raw.fraction = hdr.fraction;
records.raw.wfs = hdr.wfs;
records.raw.wfs_file = hdr.wfs;

% records.settings.wfs(1).wfs = cell(size(hdr.wfs));

% Create the first entry in the records.settings field

if param.records.file_version == 406
  check_fields = [1 4 5 6 7 8 9 10 11 12];
elseif param.records.file_version == 405
  check_fields = [1 4 5 6 7 8 9 10 11];
end

hcount = 0;
for fidx=1:size(hdr.wfs,1)
  for hidx = 1:length(hdr.wfs{fidx}.num_samp)
      hcount = hcount + 1;
      records.settings.wfs_records(hcount) = find(hdr.data_offset{fidx} >= hdr.hdr_offset{fidx}(hidx),1)+records.relative_rec_num{1}(fidx)-1;
      records.settings.wfs(hcount).wfs(1).presums = hdr.wfs{fidx,1}.presums(hidx);
      records.settings.wfs(hcount).wfs(2).presums = hdr.wfs{fidx,2}.presums(hidx);
      records.settings.wfs(hcount).wfs(1).num_sam = hdr.wfs{fidx,1}.num_samp(hidx);
      records.settings.wfs(hcount).wfs(2).num_sam = hdr.wfs{fidx,2}.num_samp(hidx);
      records.settings.wfs(hcount).wfs(1).t0 = hdr.wfs{fidx,1}.t0(hidx);
      records.settings.wfs(hcount).wfs(2).t0 = hdr.wfs{fidx,2}.t0(hidx);
      records.settings.wfs(hcount).wfs(1).bit_shifts = hdr.wfs{fidx,1}.shifts(hidx);
      records.settings.wfs(hcount).wfs(2).bit_shifts = hdr.wfs{fidx,2}.shifts(hidx);
      records.settings.wfs(hcount).wfs(1).prf = hdr.wfs{fidx,1}.prf(hidx);
      records.settings.wfs(hcount).wfs(2).prf = hdr.wfs{fidx,2}.prf(hidx);
      records.settings.wfs(hcount).wfs(1).f0 = hdr.wfs{fidx,1}.f0(hidx);
      records.settings.wfs(hcount).wfs(2).f0 = hdr.wfs{fidx,2}.f0(hidx);
      records.settings.wfs(hcount).wfs(1).f1 = hdr.wfs{fidx,1}.f1(hidx);
      records.settings.wfs(hcount).wfs(2).f1 = hdr.wfs{fidx,2}.f1(hidx);
      records.settings.wfs(hcount).wfs(1).wf_gen_clk = hdr.wfs{fidx,1}.wf_gen_clk(hidx);
      records.settings.wfs(hcount).wfs(2).wf_gen_clk = hdr.wfs{fidx,2}.wf_gen_clk(hidx);
      records.settings.wfs(hcount).wfs(1).daq_clk = hdr.wfs{fidx,1}.daq_clk(hidx);
      records.settings.wfs(hcount).wfs(2).daq_clk = hdr.wfs{fidx,2}.daq_clk(hidx);
      records.settings.wfs(hcount).wfs(1).Tpd = hdr.wfs{fidx,1}.tpd(hidx);
      records.settings.wfs(hcount).wfs(2).Tpd = hdr.wfs{fidx,2}.tpd(hidx);
      records.settings.wfs(hcount).wfs(1).tx_win = hdr.wfs{fidx,1}.tx_win(hidx);
      records.settings.wfs(hcount).wfs(2).tx_win = hdr.wfs{fidx,2}.tx_win(hidx);
      if param.records.file_version == 406
        records.settings.wfs(hcount).wfs(1).elem_slots = hdr.wfs{fidx,1}.elem_slots(hidx);
        records.settings.wfs(hcount).wfs(2).elem_slots = hdr.wfs{fidx,2}.elem_slots(hidx);
      end
      records.settings.wfs(hcount).wfs(1).blank = hdr.wfs{fidx,1}.blank(hidx);
      records.settings.wfs(hcount).wfs(1).adc_gains = hdr.wfs{fidx,1}.adc_gains(hidx,:);
      
      records.settings.wfs(hcount).wfs(2).blank = hdr.wfs{fidx,2}.blank(hidx);
      records.settings.wfs(hcount).wfs(2).adc_gains = hdr.wfs{fidx,2}.adc_gains(hidx,:);
  end
end


% records.notes = cat(2,sprintf('\nEPRI NOTES\n%s',epri_notes), ...
%   sprintf('\nCLOCK NOTES\n%s',clock_notes), ...
%   sprintf('\nRADAR TIME NOTES\n%s',radar_time_notes));
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

