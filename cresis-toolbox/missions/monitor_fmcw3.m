% script monitor_fmcw
%
% This script processes data while on the plane. It does much of the setup
% (e.g. creating a GPS file and the param structure). It then does all
% the steps for processing and posting.
%    a specific region
%
% Author: John Paden

physical_constants;
close all;
tstart = tic;
param = [];

% =======================================================================
% User Settings
% =======================================================================
param.radar_name = 'snow3'; % param.radar_name: THIS OFTEN NEEDS TO BE SET
if strcmpi(param.radar_name,'snow3')
  param.base_path = '/cresis/scratch2/yan/NRL_testflight_data/snow/';
  xml_fn = '';
  param.adc = 1; % Specify an adc to grab files from
end

% .gps_fn_dir = Optional GPS file directory (leave empty to disable)
param.gps_fn_dir = '/cresis/snfs1/dataproducts/metadata/2014_Alaska_TOnrl/';
in_fns = {fullfile(param.gps_fn_dir,'sbet_Mission1_20140308.out')};
% gps_file_type = string containing type (csv, nmea, or applanix)
gps_file_type = 'applanix';

gps_time_offset = 2;
param_day_seg = '20140308';

param_fn = ct_filename_param('snow_param_2014_Greenland_P3.xls');

% =======================================================================
% =======================================================================
% Automated Section
% =======================================================================
% =======================================================================

params = read_param_xls(param_fn,'','post');
params = params(end);

% =======================================================================
% Load data
% =======================================================================
fprintf('========================================================\n');
fprintf('Loading data (%.1f sec)\n', toc(tstart));

% Find the files that were created after this XML file
file_prefix = sprintf('%s',param.radar_name);
fns = get_filenames(param.base_path, param.radar_name, '', '.bin');
if isempty(fns)
  fprintf('  Could not find any files which match\n');
  return;
end

for fn_idx = [1 2:100:length(fns)-1 length(fns)]
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
elseif strcmpi(gps_file_type, 'applanix')
  file_type = {'applanix'};
  params = {struct('time_reference','utc', ...
    'year',str2double(param_day_seg(1:4)), ...
    'month',str2double(param_day_seg(5:6)), ...
    'day',str2double(param_day_seg(7:8)))};
  gps_source = {'applanix-field'};
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
params(1).records.frame_mode = 2;

vector_file_idxs = file_idxs;
params(1).vectors.file.start_idx = vector_file_idxs(1);
params(1).vectors.file.stop_idx = vector_file_idxs(end);
params(1).vectors.file.base_dir = param.base_path;
params(1).vectors.file.file_prefix = param.radar_name;
params(1).vectors.file.adc_folder_name = '';

params(1).get_heights.lever_arm_fh = '';

clear('param_override');
param_override = [];
param_override.sched.type = 'no scheduler';
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
params(1).post.in_path = '';
params(1).post.echo.elev_comp = 1;
params(1).post.echo.depth_axis = '[min(Surface_Depth)-3 max(Surface_Depth) + 30 ]';

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


