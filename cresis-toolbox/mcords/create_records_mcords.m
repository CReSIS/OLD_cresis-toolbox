function create_records_mcords(param, param_override)
% create_records_mcords(param, param_override)
%
% Function for creating records file for MCORDS data. This function can
% be called as a script by:
%  1. Commenting out the function line
%  2. Setting the default param structure
%  3. Uncommenting param = [] line
% This is useful for debugging.
%
% This function should be run after the GPS file has been created.
% For example, cresis-toobox/gps/missions/make_gps_2009_antarctica_DC8.m
%
% This function's output file is used by all other parts of the processing.
%
% param = struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Author: John Paden
%
% See also: check_records_mcords, create_records_mcords_task,
%   create_records_mcords_sync, create_records_mcords_post_sync

% =====================================================================
% General Setup
% =====================================================================

dbstack_info = dbstack;
if ~exist('param','var') || isempty(param) || length(dbstack_info) == 1
  % =====================================================================
  % Debug Setup
  % =====================================================================
  param = read_param_xls(ct_filename_param('rds_param_2009_Antarctica_TO.xls'),'20091222_01');
  
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

% =====================================================================
% Setup the scheduler
% =====================================================================

if strcmpi(param.sched.type,'custom_torque')
  global ctrl; % Make this global for convenience in debugging
  ctrl = torque_new_batch(param);
  fprintf('Torque batch: %s\n', ctrl.batch_dir);
  torque_compile('create_records_mcords_task.m',ctrl.sched.hidden_depend_funs,ctrl.sched.force_compile);
  
elseif ~strcmpi(param.sched.type,'no scheduler')
  fd = [get_filenames(param.path,'','','.m',struct('recursive',1)); ...
    get_filenames(param.path,'','','.mexa64',struct('recursive',1))];
  fd_override = [get_filenames(param.path_override,'','','.m',struct('recursive',1)); ...
    get_filenames(param.path_override,'','','.mexa64',struct('recursive',1))];
  
  fd = merge_filelists(fd, fd_override);

  % Remove SVN files from path
  non_svn_path = zeros(size(fd));
  for fn_idx = 1:length(fd)
    if isempty(strfind(fd{fn_idx},'.svn'))
      non_svn_path(fn_idx) = 1;
    end
  end
  fd = fd(non_svn_path);
  
  % Initialize submission ctrl structure
  global ctrl; % Make this global for convenience in debugging
  ctrl = [];
  ctrl.cmd = 'init';
  ctrl.sched = param.sched;
  ctrl.fd = fd;
  ctrl = create_task(ctrl);
  
  % Prepare submission ctrl structure for queing jobs
  ctrl.cmd = 'task';
end

% =====================================================================
% Get the files
% =====================================================================

clear wfs hdrs;
retry_fields = {};
for adc_idx = 1:length(param.records.file.adcs)
  adc = param.records.file.adcs(adc_idx);
  
  % =====================================================================
  % Get the list of files to include in this records file
  % =====================================================================
  [base_dir,adc_folder_name,filenames{adc_idx},file_idxs] = get_segment_file_list(param,adc);

  % Log file for each ADC
  param.log_fn = ct_filename_tmp(param,'','records',sprintf('log_adc%d.txt',adc));
  
  % =================================================================
  %%% Execute tasks/jobs
  fh = @create_records_mcords_task;
  arg{1} = filenames{adc_idx};
  arg{2} = param;
  arg{3} = file_idxs;
  
  if strcmp(param.sched.type,'custom_torque')
    create_task_param.conforming = true;
    create_task_param.notes = sprintf('Getting files for adc %d (%d of %d) (%s)\n', ...
      adc, adc_idx, length(param.records.file.adcs), datestr(now,'HH:MM:SS'));
    ctrl = torque_create_task(ctrl,fh,3,arg,create_task_param);
    
  elseif ~strcmp(param.sched.type,'no scheduler')
    [ctrl,job_id,task_id] = create_task(ctrl,fh,3,arg);
    fprintf('Getting files for adc %d (%d of %d) (%s)\n', ...
      adc, adc_idx, length(param.records.file.adcs), datestr(now,'HH:MM:SS'));    
    retry_fields{job_id,task_id}.adc = adc;
    retry_fields{job_id,task_id}.adc_idx = adc_idx;
    retry_fields{job_id,task_id}.length = length(param.records.file.adcs);
    retry_fields{job_id,task_id}.arg = arg;
  
  else
    fprintf('Getting files for adc %d (%d of %d) (%s)\n', ...
      adc, adc_idx, length(param.records.file.adcs), datestr(now,'HH:MM:SS'));    
    [success tmp_hdr wfs{adc_idx}] = fh(arg{1},arg{2},arg{3});
    hdrs.ver{adc_idx} = tmp_hdr.ver;
    hdrs.seconds{adc_idx} = tmp_hdr.seconds;
    hdrs.fractions{adc_idx} = tmp_hdr.fractions;
    hdrs.epri{adc_idx} = tmp_hdr.epri;
    hdrs.offset{adc_idx} = tmp_hdr.offset;
    hdrs.file_idx{adc_idx} = tmp_hdr.file_idx;
  end
  
end

% =======================================================================
% Wait for jobs to complete if a scheduler was used
% =======================================================================
if strcmpi(param.sched.type,'custom_torque')
  % Wait until all submitted jobs to complete
  ctrl = torque_rerun(ctrl);
  if ~all(ctrl.error_mask == 0)
    if ctrl.sched.stop_on_fail
      torque_cleanup(ctrl);
      error('Not all jobs completed, but out of retries (%s)', datestr(now));
    else
      warning('Not all jobs completed, but out of retries (%s)', datestr(now));
      keyboard;
    end
  else
    fprintf('Jobs completed (%s)\n\n', datestr(now));
  end
  for adc_idx = 1:length(param.records.file.adcs)
    out_fn = fullfile(ctrl.out_path_dir,sprintf('out_%d.mat',adc_idx));
    out = load(out_fn);
    hdrs.ver{adc_idx} = out.argsout{2}.ver;
    hdrs.seconds{adc_idx} = out.argsout{2}.seconds;
    hdrs.fractions{adc_idx} = out.argsout{2}.fractions;
    hdrs.epri{adc_idx} = out.argsout{2}.epri;
    hdrs.offset{adc_idx} = out.argsout{2}.offset;
    hdrs.file_idx{adc_idx} = out.argsout{2}.file_idx;
    wfs{adc_idx} = out.argsout{3};
  end
  torque_cleanup(ctrl);
  
elseif ~strcmpi(param.sched.type,'no scheduler')
  ctrl.cmd = 'done';
  ctrl = create_task(ctrl);
  fprintf('Jobs completed (%s)\n\n', datestr(now,'HH:MM:SS'));
  
  % No errors, so get outputs
  for adc_idx = 1:length(param.records.file.adcs)
    adc_idx
    wfs{adc_idx} = ctrl.jobs{job_idxs(adc_idx)}.args_out{task_idxs(adc_idx),3};
    hdrs.ver{adc_idx} = ctrl.jobs{job_idxs(adc_idx)}.args_out{task_idxs(adc_idx),2}.ver;
    hdrs.seconds{adc_idx} = ctrl.jobs{job_idxs(adc_idx)}.args_out{task_idxs(adc_idx),2}.seconds;
    hdrs.fractions{adc_idx} = ctrl.jobs{job_idxs(adc_idx)}.args_out{task_idxs(adc_idx),2}.fractions;
    hdrs.epri{adc_idx} = ctrl.jobs{job_idxs(adc_idx)}.args_out{task_idxs(adc_idx),2}.epri;
    hdrs.offset{adc_idx} = ctrl.jobs{job_idxs(adc_idx)}.args_out{task_idxs(adc_idx),2}.offset;
    hdrs.file_idx{adc_idx} = ctrl.jobs{job_idxs(adc_idx)}.args_out{task_idxs(adc_idx),2}.file_idx;
  end
end

% Count the presums
num_presum = 0;
for wf = 1:length(wfs{1}.presums)
  num_presum = num_presum + wfs{1}.presums(wf);
end

%% Save workspace in case there is a failure
create_records_save_workspace;

%% Correct time, sync GPS data, and save records
create_records_mcords_sync;

return;
