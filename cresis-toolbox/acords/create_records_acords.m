function create_records_acords(param,param_override)
% create_records_acords(param,param_override)
%
% Corrects jumps in utc time, typically one second jumps. This script
% obtains fmcw headers from data from data drives indicated in the param
% spreadsheets. After loading the files using basic_load.m the header files
% are saved in the support directories for the specific radar in .mat form
% for quicker access.
%
% Can be run as a function by passing in the param argument
% or a script (by setting the default value of param).
%
% Output file contains:
% hdr: structure with the following fields
%    .utc_time_sod: vector of corrected utc times
%    .seconds: vector
%    .fraction: vector
%
% Inputs:
% param = struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Authors: Aric Beaver, John Paden, Logan Smith
%
% See also: master

% =====================================================================
% General Setup
% =====================================================================

dbstack_info = dbstack;
if ~exist('param','var') || isempty(param) || length(dbstack_info) == 1
  % =====================================================================
  % Debug Setup
  % =====================================================================
  param = read_param_xls(ct_filename_param('rds_param_2005_Greenland_TO.xls'),'20050505_01');
  
  clear('param_override');
  param_override.sched.type = 'no scheduler';
  param_override.sched.rerun_only = true;
  
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
  
elseif ~isstruct(param)
  % Functional form
  param();
end
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

% =====================================================================
%% Check input parameters
% =====================================================================

if ~isfield(param.records,'force_all') || isempty(param.records.force_all)
  param.records.force_all = false;
end

if ~isfield(param.records,'file_version') || isempty(param.records.file_version)
  param.records.file_version = 406;
end

if ~isfield(param.records,'debug_level') || isempty(param.records.debug_level)
  param.records.debug_level = 1;
end

param.sched.type = 'no scheduler';

% =====================================================================
%% Get the list of raw files to load in
% =====================================================================

%%%%%%%%%%%%%%%% Add acords case to get_segment_file_list
[base_dir,adc_folder_name,fns,file_idxs] = get_segment_file_list(param);

% =====================================================================
%% Load headers from radar data files
% =====================================================================

% Load header from radar data file
fprintf('Loading raw files %i to %i\n',file_idxs([1 end]));

clear hdr
hdr = [];
hdr.seconds = zeros([0 0],'double');
hdr.fraction = zeros([0 0],'uint32');
% if param.records.file_version ~= 101
%   hdr.epri = zeros([0 0],'uint32');
% end
hdr.offset = zeros([0 0],'uint32');
hdr.wfs = [];
hdr.wfs_file = [];
hdr.hdr_offset = cell(1,length(file_idxs));
hdr.data_offset = cell(1,length(file_idxs));
% if param.records.file_version == 2
%   hdr.nyquist_zone = zeros([0 0],'uint8');
%   hdr.loopback_mode = zeros([0 0],'uint8');
% end
records = [];
% records.relative_rec_num = This variable contains the first record
% number of each file. After the loop runs there will always be one
% to many elements (161 files will mean 162 elements in the array)
% and the last entry is the record number that would have been next
% so that length(hdr.utc_time_sod) = records.relative_rec_num(end)-1
records.relative_rec_num{1} = 1;

%init_EPRI_estimate = create_records_epri_estimate(param,file_idxs,fns);

for file_idx = 1:length(file_idxs)
  file_num = file_idxs(file_idx);
  fn = fns{file_num};
  
  %% Prepare arguments based on param.radar_name
  first_byte = 0;
%   arg{1} = fns{file_num};
%   arg{2} = struct('clk',param.radar.fs,'utc_time_halved',param.vectors.gps.utc_time_halved, ...
%     'first_byte',first_byte, 'file_version', param.records.file_version, ...
%     'records',struct('en',1,'epri',init_EPRI_estimate,'force_all',param.records.force_all));
%   fh = @basic_load_acords;
  
  %% Create/run task for each file
  fprintf('  %i/%i %s (%s)\n', ...
    file_idx,length(file_idxs), fns{file_num}, datestr(now,'HH:MM:SS'));
  
  [~,fn_name,ext] = fileparts(fn);
  
  if strcmp(fn_name,'nov18')
    fn_name = 'nov18_04';
  end
  
  %     if isfield(param.records,'tmp_fn_uses_adc_folder_name') && param.records.tmp_fn_uses_adc_folder_name
  tmp_hdr_fn = ct_filename_ct_tmp(rmfield(param,'day_seg'),'','headers', ...
    fullfile(adc_folder_name, [fn_name ext '.mat']));
  %     else
  %       tmp_hdr_fn = ct_filename_ct_tmp(rmfield(param,'day_seg'),'','headers',[fn_name ext '.mat']);
  %     end
  if exist(tmp_hdr_fn,'file')
    clear hdr_tmp
    hdr_tmp = load(tmp_hdr_fn);
  else
    error('%s does not exist!',tmp_hdr_fn)
  end
  % Concatenate hdr_tmp fields
  hdr.seconds = cat(2,hdr.seconds,hdr_tmp.seconds);
  %     hdr.fraction = cat(2,hdr.fraction,hdr_tmp.fraction);
  %     if length(hdr.seconds) ~= length(hdr.fraction)
  %       keyboard
  %     end
  hdr.offset = cat(2,hdr.offset,hdr_tmp.offset);
  hdr.wfs = cat(1,hdr.wfs,hdr_tmp.wfs);  
  hdr.hdr_offset{file_idx} = hdr_tmp.hoffset;
  hdr.data_offset{file_idx} = hdr_tmp.offset;
  
  
  %% Create records and file numbers
  records.relative_rec_num{1}(file_idx+1) = length(hdr_tmp.seconds)+records.relative_rec_num{1}(file_idx);
  [fn_dir fn_name fn_ext] = fileparts(fn);
  records.relative_filename{1}{file_idx} = [fn_name fn_ext];
  hdr.wfs_file = cat(2,hdr.wfs_file,file_idx*ones(1,length(hdr_tmp.hdr)));
  
end
records.relative_rec_num{1} = records.relative_rec_num{1}(1:end-1);

hdr.fraction = zeros(size(hdr.seconds));

%% Save workspace in case there is a failure
create_records_save_workspace;

%% Correct time, sync GPS data, and save records
create_records_acords_sync;

return;

