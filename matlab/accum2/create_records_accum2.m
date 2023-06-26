function create_records_accum2(param, param_override)
% create_records_accum2(param, param_override)
%
% Function for creating records file for accum2 data. This function can
% be called as a script by:
%  1. Commenting out the function line
%  2. Setting the default param structure
%  3. Uncommenting param = [] line
% This is useful for debugging.
%
% This function should be run after the GPS file has been created.
% For example, cresis-toobox/gps/missions/make_gps_2009_antarctica_DC8.m
%
% This function's output is used by create_frames.m.
%
% Author: John Paden
%
% See also: check_records_mcrds, create_records_mcrds_task,create_records_mcrds_post_sync

% =====================================================================
% General Setup
% =====================================================================

dbstack_info = dbstack;
if ~exist('param','var') || isempty(param) || length(dbstack_info) == 1
  % =====================================================================
  % Debug Setup
  % =====================================================================
  param = read_param_xls(ct_filename_param('accum_param_2013_Antarctica_Ground.xls'),'20131219_02');
  
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
% Prep work
% =====================================================================

% Get WGS84 ellipsoid parameters
physical_constants

radar_name = param.radar_name;
ver = 1;

% =====================================================================
% Get the files
% =====================================================================

clear wfs hdrs job_idxs task_idxs;
  
fprintf('Getting files (%s)\n', datestr(now,'HH:MM:SS'));
  
% =====================================================================
% Get the list of files to include in this records file
% =====================================================================

[base_dir,adc_folder_name,fns,file_idxs] = get_segment_file_list(param);

% Read in the header information from all of the files
records = create_records_accum2_task(fns,file_idxs);
if strcmpi(param.season_name,'2013_Antarctica_P3') & strcmpi(param.day_seg,'20131120_01') 
  records.num_coh_ave(:) = 102;
  records.wfs{1}.presums = 102;
end

% Count the presums in one EPRI (assumes number of waveforms
% and number of presums does not change)
% num_presum = records.wfs{1}.num_coh_ave(1) * length(records.wfs);
num_presum = records.num_coh_ave(1) * length(records.wfs);

% ===================================================================
% ===================================================================
% Synchronize radar data to GPS data
% ===================================================================
% ===================================================================
if isfield(param,'tmp_path') && ~isempty(param.tmp_path)
  fn = ct_filename_tmp(param,'','records','workspace');
  fprintf('Saving workspace %s\n', fn);
  try
    fn_dir = fileparts(fn);
    if ~exist(fn_dir,'dir')
      mkdir(fn_dir);
    end
    save(fn);
  catch
    fprintf('  Saving workspace failed\n');
  end
end

create_records_accum2_post;

return;
