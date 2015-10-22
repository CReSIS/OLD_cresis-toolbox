% script monitor_mcords4
%
% This script processes data while on the plane. It does much of the setup
% (e.g. creating a GPS file and the param structure). It then does all
% the steps for processing and posting.
%    a specific region
%
% Author: John Paden

error('Use monitor mcords3 since it should work for mcords4 too. Once tested, remove this function.')
physical_constants;
close all;
tstart = tic;
param = [];

% =======================================================================
% User Settings
% =======================================================================
param.radar_name = 'mcords4'; % param.radar_name: THIS OFTEN NEEDS TO BE SET
if strcmpi(param.radar_name,'mcords4')
  param.base_path = '/basler/fp1/mcords4/';
  xml_fn = '';
  param.adc = 1; % Specify an adc to grab files from
end

% .gps_fn_dir = Optional GPS file directory (leave empty to disable)
param.gps_fn_dir = '/basler/fp1/radar_gps/';
in_fns = {fullfile(param.gps_fn_dir,'gps_20140107.csv')};

gps_time_offset = 2-86400;
param_day_seg = '20140108_01';

param_fn = '/basler/scratch2/scripts/params-cr1/rds_param_2013_Antarctica_Basler.xls';

% =======================================================================
% =======================================================================
% Automated Section
% =======================================================================
% =======================================================================

% =======================================================================
% Load data
% =======================================================================
fprintf('========================================================\n');
fprintf('Loading data (%.1f sec)\n', toc(tstart));


% First, get the last XML file made
if isempty(xml_fn)
  xml_fns = get_filenames(param.base_path,'','','.xml');
  for xml_idx = 1:length(xml_fns)
    xml_fn = xml_fns{xml_idx};
    fprintf('=================== xml index %d ================\n', xml_idx);
    fprintf('%s\n', xml_fn);
    [settings,settings_enc] = read_cresis_xml(xml_fn);
    fprintf(' # of waveforms: %d\n', settings.DDS_Setup.Wave);
    settings.DDS_Setup.Waveforms(1)
  end
  xml_idx = input('Input xml index: ');
  xml_fn = xml_fns{xml_idx};
end

[settings,settings_enc] = read_cresis_xml(xml_fn);

xml_fns = get_filenames(param.base_path,'','','.xml');
for xml_idx = 1:length(xml_fns)
  xml_fn = xml_fns{xml_idx};
  [tmp xml_fn_name] = fileparts(xml_fn);
  time_stamp_idx = find(xml_fn_name == '_',1) + 1;
  year = str2double(xml_fn_name(time_stamp_idx + (0:3)));
  month = str2double(xml_fn_name(time_stamp_idx + 4 + (0:1)));
  day = str2double(xml_fn_name(time_stamp_idx + 6 + (0:1)));
  hour = str2double(xml_fn_name(time_stamp_idx + 9 + (0:1)));
  min = str2double(xml_fn_name(time_stamp_idx + 11 + (0:1)));
  sec = str2double(xml_fn_name(time_stamp_idx + 13 + (0:1)));
  %fprintf('  %s: %d\n', new_settings.fn, length(new_settings.Configuration.Waveforms));
  xml_date(xml_idx) = datenum(year,month,day,hour,min,sec);
end

xml_idx = find(strcmp(settings.fn, xml_fns));
date_begin = settings.datenum;
if xml_idx == length(xml_fns)
  date_end = inf;
else
  date_end = xml_date(xml_idx+1);
end

% Find the files that were created after this XML file
param.board = adc_to_board('mcords4',param.adc);
file_prefix = sprintf('mcords4_%02d_',param.board);
fns = get_filenames(fullfile(param.base_path,sprintf('chan%d',param.board)), 'mcords4', '', '.bin');
fn_dir = fullfile(param.base_path, sprintf('chan%d',param.board));
if isempty(fns)
  fprintf('  Could not find any files which match\n');
  return;
end

fn_mask = logical(size(fns));
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

gps_path = fullfile(support_path,'gps','2013_Antarctica_Basler');
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

file_type = {'CSV'};
params = {struct('time_reference','utc')};
gps_source = {'csv-field'};
sync_flag = {0};
make_gps;

% =========================================================================
% =========================================================================
% 2. Use param sheet to get all parameters, except over ride vector files

params = read_param_xls(param_fn,'','post');
params = params(end);
params(1).day_seg = param_day_seg;
params(1).vectors.gps.time_offset = gps_time_offset;
params(1).records.frame_mode = 2;

vector_file_idxs = find(fn_mask);
vector_file_idxs = vector_file_idxs(file_idxs);
params(1).vectors.file.start_idx = vector_file_idxs(1);
params(1).vectors.file.stop_idx = vector_file_idxs(end);
params(1).vectors.file.base_dir = param.base_path;

params(1).radar.prf = settings.DDS_Setup.PRF;
params(1).prf = settings.DDS_Setup.PRF;
param_wf = params(1).radar.wfs(1);
params(1).radar.wfs = param_wf;
for wf = 1:length(settings.DDS_Setup.Waveforms)
  params(1).radar.wfs(wf) = param_wf;
  params(1).radar.wfs(wf).Tpd = double(settings.DDS_Setup.Waveforms(wf).Len_Mult) ...
    * settings.DDS_Setup.Base_Len;
  params(1).radar.wfs(wf).f0 = settings.DDS_Setup.Waveforms(wf).Start_Freq(1);
  params(1).radar.wfs(wf).f1 = settings.DDS_Setup.Waveforms(wf).Stop_Freq(1);
  params(1).radar.wfs(wf).tukey = settings.DDS_Setup.RAM_Taper;
  params(1).radar.wfs(wf).tx_weights = double(settings.DDS_Setup.Ram_Amplitude)/60000*250;
  params(1).radar.wfs(wf).adc_gains = 10.^((52-double(settings.DDS_Setup.Waveforms(wf).Attenuator_1 + settings.DDS_Setup.Waveforms(wf).Attenuator_2))/20);
end
Tpd = cell2mat({params(1).radar.wfs(:).Tpd});
params(1).get_heights.qlook.img_comb = [Tpd(3:2:end); Tpd(1:2:end-2)];
params(1).get_heights.qlook.img_comb = reshape(params(1).get_heights.qlook.img_comb,[1 numel(params(1).get_heights.qlook.img_comb)]);

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
params(1).post.echo.depth_axis = '[min(Surface_Elev)-1500 max(Surface_Elev) + 100 ]';
params(1).post.img_dpi = 100;

wfs = 1:2:length(settings.DDS_Setup.Waveforms);
Tpd= double(cell2mat({settings.DDS_Setup.Waveforms(wfs).Len_Mult}))*settings.DDS_Setup.Base_Len;
[~,sort_idx] = sort(Tpd);
    
img_comb = [];
if length(settings.DDS_Setup.Waveforms) > 2
  wf_idx = 2;
  Tpd_prev = double(settings.DDS_Setup.Waveforms(wfs(wf_idx-1)).Len_Mult)*settings.DDS_Setup.Base_Len;
  Tpd = double(settings.DDS_Setup.Waveforms(wfs(wf_idx)).Len_Mult)*settings.DDS_Setup.Base_Len;
  img_comb = cat(2,img_comb,[Tpd Tpd_prev]);
  for wf_idx = 3:length(wfs)
    Tpd_prev = double(settings.DDS_Setup.Waveforms(wfs(wf_idx-1)).Len_Mult)*settings.DDS_Setup.Base_Len;
    Tpd = double(settings.DDS_Setup.Waveforms(wfs(wf_idx)).Len_Mult)*settings.DDS_Setup.Base_Len;
    img_comb = cat(2,img_comb,[Tpd Tpd_prev]);
  end
end

imgs = {};

for wf_idx = 1:length(wfs)
  wf = wfs(wf_idx);
  imgs{wf_idx} = [];
  for adc = 1:length(settings.DDS_Setup.Waveforms(1).Start_Freq)
    imgs{wf_idx} = cat(1,imgs{wf_idx},[-j*wf adc]);
  end
end


if 0
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


