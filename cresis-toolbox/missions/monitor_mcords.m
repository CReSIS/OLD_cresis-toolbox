% script monitor_mcords3
%
% This script processes data while on the plane. It does much of the setup
% (e.g. creating a GPS file and the param structure). It then does all
% the steps for processing and posting.
%
% Author: John Paden

physical_constants;
close all;
tstart = tic;
param = [];

% =======================================================================
% User Settings
% =======================================================================
param.radar_name = 'mcords3'; % param.radar_name: THIS OFTEN NEEDS TO BE SET
if strcmpi(param.radar_name,'mcords3')
  param.base_path = '/process/mcords/';
  xml_fn = '';
  param.adc = 1; % Specify an adc to grab files from
  iq_mode = 0;
  
  % xml_version = There are three versions of XML files
  %   mcords3: pre-2014 Greenland P3 is 1
  %   mcords3: post-2014 Greenland P3 is 3
  %   mcords4: 2
  xml_version = 3;
  
elseif strcmpi(param.radar_name,'mcords4')
  param.base_path = '/basler/fp1/mcords4/';
  xml_fn = '';
  param.adc = 1; % Specify an adc to grab files from
  iq_mode = 1;

  % xml_version = There are three versions of XML files
  %   mcords3: pre-2014 Greenland P3 is 1
  %   mcords3: post-2014 Greenland P3 is 3
  %   mcords4: 2
  xml_version = 2;
end

% .gps_fn_dir = Optional GPS file directory (leave empty to disable)
param.gps_fn_dir = '/scratch/metadata/2014_Greenland_P3/';
in_fns = {fullfile(param.gps_fn_dir,'GPS_030614_195145.txt')};
%gps_file_type = 'csv';
gps_file_type = 'nmea';

gps_time_offset = 1;
param_day_seg = '20140307_01';

param_fn = ct_filename_param('rds_param_2014_Greenland_P3.xls');

% =======================================================================
% =======================================================================
% Automated Section
% =======================================================================
% =======================================================================

params = read_param_xls(param_fn,'','post');
params = params(end);
if strcmpi(param.radar_name,'mcords3')
  if xml_version == 1
    config_var = 'Configuration';
    config_var_enc = 'Configuration';
    prf_var = 'PRF_Hz';
    ram_var = 'RAM Taper';
    ram_var_enc = 'RAMZ20Taper';
    xml_file_prefix = 'DDS';
    phase_var = 'Phase_Offset_deg';
    phase_var_enc = 'PhaseZ20OffsetZ20Z28degZ29';
    wave_var_enc = 'Z23Waveforms';
    ttl_start_var_enc = 'TTLZ20StartZ20Z28TTLZ30Z2DZ3ETTLZ37Z29';
    ttl_length_var_enc = 'TTLZ20LengthZ20Z28TTLZ30Z2DZ3ETTLZ37Z29';
    TTL_prog_delay = 650;
  elseif xml_version == 3
    config_var = 'DDS_Setup';
    config_var_enc = 'DDSZ20Setup';
    prf_var = 'PRF';
    ram_var = 'Ram_Amplitude';
    ram_var_enc = 'RamZ20Amplitude';
    xml_file_prefix = 'radar';
    phase_var = 'Phase_Offset';
    phase_var_enc = 'PhaseZ20Offset';
    wave_var_enc = 'Z23Wave';
    ttl_start_var_enc = 'TTLZ20Start';
    ttl_length_var_enc = 'TTLZ20Length';
    TTL_prog_delay = 650;
  else
    error('Unsupported version');
  end
  fs = 1e9/9;
  fs_sync = 1e9/18;
elseif strcmpi(param.radar_name,'mcords4')
  if xml_version == 2
    config_var = 'DDS_Setup';
    config_var_enc = 'DDSZ20Setup';
    prf_var = 'PRF';
    ram_var = 'Ram_Amplitude';
    ram_var_enc = 'RamZ20Amplitude';
    xml_file_prefix = 'mcords4';
    phase_var = 'Phase_Offset';
    phase_var_enc = 'PhaseZ20Offset';
    wave_var_enc = 'Z23Wave';
    ttl_start_var_enc = 'TTLZ20Start';
    ttl_length_var_enc = 'TTLZ20Length';
    TTL_prog_delay = 650;
    fs = 1e9/2;
    fs_sync = 1e9/8;
  else
    error('Unsupported version');
  end
end

% =======================================================================
% Load data
% =======================================================================
fprintf('========================================================\n');
fprintf('Loading data (%.1f sec)\n', toc(tstart));

% First, get the last XML file made
if isempty(xml_fn)
  [settings,settings_enc] = read_ni_xml_directory(param.base_path,xml_file_prefix,false);
  for xml_idx = 1:length(settings)
    xml_fn = settings(xml_idx).fn;
    fprintf('=================== xml index %d ================\n', xml_idx);
    fprintf('%s\n', xml_fn);
    fprintf(' # of waveforms: %d\n', length(settings(xml_idx).(config_var).Waveforms));
    settings(xml_idx).(config_var).Waveforms(1)
  end
  xml_idx = input('Input xml index: ');
  xml_fn = settings(xml_idx).fn;
end
  
% xml_idx = find(strcmp(settings.fn, xml_fns));
date_begin = settings(xml_idx).datenum;
if xml_idx == length(settings)
  date_end = inf;
else
  date_end = settings(xml_idx+1).datenum;
end
settings = settings(xml_idx);
settings_enc = settings_enc(xml_idx);
for wf = 1:length(settings.(config_var).Waveforms)
  if ~isfield(settings.(config_var).Waveforms(wf), 'Delay')
    settings.(config_var).Waveforms(wf).Delay = zeros(size(settings.(config_var).Waveforms(wf).(phase_var)));
  end
end

% Find the files that were created after this XML file
param.board = adc_to_board(param.radar_name,param.adc);
file_prefix = sprintf('%s_%02d_',param.radar_name,param.board);
if strcmpi(param.radar_name,'mcords3')
  fns = get_filenames(fullfile(param.base_path,sprintf('board%d',param.board)), 'mcords3', '', '.bin');
  fn_dir = fullfile(param.base_path, sprintf('board%d',param.board));
elseif strcmpi(param.radar_name,'mcords4')
  fns = get_filenames(fullfile(param.base_path,sprintf('chan%d',param.board)), 'mcords4', '', '.bin');
  fn_dir = fullfile(param.base_path, sprintf('chan%d',param.board));
end
if isempty(fns)
  fprintf('  Could not find any files which match\n');
  return;
end

fn_mask = logical(zeros(size(fns)));
for fn_idx = 1:length(fns)
  fn = fns{fn_idx};
  
  fname = fname_info_mcords2(fn);
  fn_date(fn_idx) = fname.datenum;
  
  if fn_date(fn_idx) >= date_begin && fn_date(fn_idx) <= date_end
    fn_mask(fn_idx) = 1;
  end
end

fns = fns(fn_mask);
for fn_idx = 1:length(fns)
  fprintf('%d: %s\n', fn_idx, fns{fn_idx});
end
file_idxs = input('Input file indexes: ');
fns = fns(file_idxs);

% =========================================================================
% =========================================================================
%% 1. Make GPS file
tic;
global gRadar;

support_path = '';
data_support_path = '';

if isempty(support_path)
  support_path = gRadar.support_path;
end

out_fns = {sprintf('gps_%s.mat',param_day_seg(1:8))};

gps_path = fullfile(support_path,'gps',params.season_name);
if ~exist(gps_path,'dir')
  fprintf('Making directory %s\n', gps_path);
  fprintf('  Press a key to proceed\n');
  pause;
  mkdir(gps_path);
end

if isempty(data_support_path)
  data_support_path = gRadar.data_support_path;
end

debug_level = 1;

if strcmpi(gps_file_type, 'csv')
  file_type = {'CSV'};
  params = {struct('time_reference','utc')};
  gps_source = {'csv-field'};
  sync_flag = {0};
elseif strncmpi(gps_file_type, 'nmea',4)
  file_type = {'NMEA'};
  [~,in_fns_name] = fileparts(in_fns{1});
  if strncmpi(in_fns_name,'accum',5)
    nmea_format = 3;
  else
    nmea_format = 1;
  end
  params = {struct('time_reference','utc', ...
    'year',str2double(param_day_seg(1:4)), ...
    'month',str2double(param_day_seg(5:6)), ...
    'day',str2double(param_day_seg(7:8)), ...
    'format',nmea_format, 'nmea_tag', '$GPGGA')};
  gps_source = {'nmea-field'};
  sync_flag = {0};
end
make_gps;

% =========================================================================
% =========================================================================
% 2. Use param sheet to get all parameters, except over ride vector files

params = read_param_xls(param_fn,'','post');
params = params(end);
params(1).day_seg = param_day_seg;
params(1).vectors.gps.time_offset = gps_time_offset;
params(1).records.frame_mode = bin2dec('101');

vector_file_idxs = find(fn_mask);
vector_file_idxs = vector_file_idxs(file_idxs);
params(1).vectors.file.start_idx = vector_file_idxs(1);
params(1).vectors.file.stop_idx = vector_file_idxs(end);
params(1).vectors.file.base_dir = param.base_path;
if strcmpi(param.radar_name,'mcords3')
  params(1).vectors.file.adc_folder_name = 'board%b';
elseif strcmpi(param.radar_name,'mcords4')
  params(1).vectors.file.adc_folder_name = 'chan%d';
end

params(1).radar.prf = settings.(config_var).(prf_var);
params(1).prf = settings.(config_var).(prf_var);
param_wf = params(1).radar.wfs(1);
params(1).radar.wfs = param_wf;
for wf = 1:length(settings.(config_var).Waveforms)
  params(1).radar.wfs(wf) = param_wf;
  params(1).radar.wfs(wf).Tpd = double(settings.(config_var).Waveforms(wf).Len_Mult) ...
    * settings.(config_var).Base_Len;
  params(1).radar.wfs(wf).f0 = settings.(config_var).Waveforms(wf).Start_Freq(1);
  params(1).radar.wfs(wf).f1 = settings.(config_var).Waveforms(wf).Stop_Freq(1);
  params(1).radar.wfs(wf).tukey = settings.(config_var).RAM_Taper;
  params(1).radar.wfs(wf).tx_weights = double(settings.(config_var).(ram_var))/60000*250;
  params(1).radar.wfs(wf).adc_gains = 10.^((52-double(settings.(config_var).Waveforms(wf).Attenuator_1 + settings.(config_var).Waveforms(wf).Attenuator_2))/20);
end
Tpd = cell2mat({params(1).radar.wfs(:).Tpd});
params(1).get_heights.qlook.img_comb = [Tpd(3:2:end); Tpd(1:2:end-2)];
params(1).get_heights.qlook.img_comb = reshape(params(1).get_heights.qlook.img_comb,[1 numel(params(1).get_heights.qlook.img_comb)]);
params(1).get_heights.lever_arm_fh = '';

clear('param_override');
param_override = [];
% param_override.sched.type = 'no scheduler';
param_override.sched.cluster_size = inf;
param_override.sched.rerun_only = false;
param_override.sched.submit_arguments    = '-l nodes=1:ppn=1,walltime=60:00';
param_override.sched.stop_on_fail = false;

% params(1).cmd.create_vectors = 0;
% params(1).cmd.create_records = 0;
% params(1).cmd.create_frames = 0;
params(1).cmd.create_vectors = 1;
params(1).cmd.create_records = 1;
params(1).cmd.create_frames = 1;
params(1).cmd.get_heights = 1;
params(1).cmd.csarp = 0;
params(1).cmd.combine_wf_chan = 0;
params(1).cmd.generic = 1;
params(1).post.echo.elev_comp = 3;
params(1).post.echo.depth_axis = '[min(Surface_Elev)-3500 max(Surface_Elev) + 100 ]';
params(1).post.img_dpi = 100;

if abs(iq_mode)
  wfs = 1:2:length(settings.(config_var).Waveforms);
else
  wfs = 1:length(settings.(config_var).Waveforms);
end
Tpd= double(cell2mat({settings.(config_var).Waveforms(wfs).Len_Mult}))*settings.(config_var).Base_Len;
[~,sort_idx] = sort(Tpd);
    
img_comb = [];
if length(settings.(config_var).Waveforms) > 1+abs(iq_mode)
  wf_idx = 2;
  Tpd_prev = double(settings.(config_var).Waveforms(wfs(wf_idx-1)).Len_Mult)*settings.(config_var).Base_Len;
  Tpd = double(settings.(config_var).Waveforms(wfs(wf_idx)).Len_Mult)*settings.(config_var).Base_Len;
  img_comb = cat(2,img_comb,[Tpd Tpd_prev]);
  for wf_idx = 3:length(wfs)
    Tpd_prev = double(settings.(config_var).Waveforms(wfs(wf_idx-1)).Len_Mult)*settings.(config_var).Base_Len;
    Tpd = double(settings.(config_var).Waveforms(wfs(wf_idx)).Len_Mult)*settings.(config_var).Base_Len;
    img_comb = cat(2,img_comb,[Tpd Tpd_prev]);
  end
end

imgs = {};

for wf_idx = 1:length(wfs)
  wf = wfs(wf_idx);
  imgs{wf_idx} = [];
  for adc = 1:length(settings.(config_var).Waveforms(1).Start_Freq)
    if abs(iq_mode)
      imgs{wf_idx} = cat(1,imgs{wf_idx},[iq_mode*wf adc]);
    else
      imgs{wf_idx} = cat(1,imgs{wf_idx},[wf adc]);
    end
  end
end


if 1
  params(1).get_heights.imgs = imgs;
  params(1).get_heights.qlook.img_comb = img_comb;
  param_override.sched.rerun_only = false;
else
  params(1).get_heights.imgs = {[-1j 1; -1j 2; -1j 3; -1j 4; -1j 5; -1j 6; -1j 7; -1j 8]}
  params(1).get_heights.qlook.img_comb = [];
  param_override.sched.rerun_only = false;
end

% =========================================================================
% =========================================================================
% 3. Run everything
master(params,param_override);

% =========================================================================
% =========================================================================
% 3. Post results
params(1).post.data_dirs = {'qlook'};
params(1).post.layer_dir = '';
params(1).post.layers_en = 0;
params(1).post.data_en = 0;
params(1).post.csv_en = 0;
params(1).post.concat_en = 0;
create_posting(params,param_override);




return


