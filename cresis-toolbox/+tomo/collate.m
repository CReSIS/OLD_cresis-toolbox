function collate(param, param_override)
% tomo.collate(param, param_override)
%
% Usually this function is called from tomo.run_collate.
% Calls tomo.collate_task.
%
% Inputs:
%   param = struct with processing parameters
%   param_override = parameters in this struct will override parameters
%     in param.
%
% See also: tomo.run_collate, tomo.collate, tomo_collate_task,
%   tomo.fuse_images, tomo.add_icemask_surfacedem, tomo.create_surfData,
%
% Author: John Paden, Jordan Sprick, and Mingze Xu

param = merge_structs(param,param_override);

%% Determine which frames we will operate on
% Load frames file
if ~isfield(param.records,'frames_fn')
  param.records.frames_fn = '';
end
load(ct_filename_support(param,param.records.frames_fn,'frames'));

if isempty(param.cmd.frms)
  param.cmd.frms = 1:length(frames.frame_idxs);
end
% Remove frames that do not exist from param.cmd.frms list
[valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
% valid_frms = ones(1,length(param.cmd.frms));
if length(valid_frms) ~= length(param.cmd.frms)
  bad_mask = ones(size(param.cmd.frms));
  bad_mask(keep_idxs) = 0;
  warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
    param.cmd.frms(find(bad_mask,1)));
  param.cmd.frms = valid_frms;
end

%% Compile C++ functions
if 0
  % If you get a C++11 option error, you may be using pre-G++ 4.7. You can
  % check the g++ version with system('g++ --version');
  % To fix this, add -v option to mex function and look for a line like this:
  %   Options file: ~/.matlab/R2015b/mex_C++_glnxa64.xml
  % Replace -std=c++11 with -std=c++0x (should occur in two places)
  % Reference: http://stackoverflow.com/questions/14674597/cc1plus-error-unrecognized-command-line-option-std-c11-with-g
  mex -largeArrayDims fuse.cpp
  mex -largeArrayDims train_model.cpp
  mex -largeArrayDims detect.cpp
  mex -largeArrayDims extract.cpp
  mex -largeArrayDims refine.cpp
end

%% Initialize Torque setup
if strcmpi(param.sched.type,'custom_torque')
  global ctrl; % Make this global for convenience in debugging
  ctrl = torque_new_batch(param);
  fprintf('Torque batch: %s\n', ctrl.batch_dir);
  torque_compile('collate_task.m',ctrl.sched.hidden_depend_funs,ctrl.sched.force_compile);
end

%% Create Tasks
task_param = param;
fh = @collate_task;
for frm_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frm_idx);
  task_param.proc.frm = frm;
  arg = {task_param};
  
  if strcmp(param.sched.type,'custom_torque')
    create_task_param.conforming = true;
    create_task_param.notes = sprintf('%s_%03d (%d of %d)', ...
      param.day_seg, frm, frm_idx, length(param.cmd.frms));
    ctrl = torque_create_task(ctrl,fh,1,arg,create_task_param);
  else
    fprintf('%s_%03d (%d of %d)\n', ...
      param.day_seg, frm, frm_idx, length(param.cmd.frms));
    fh(arg{:});
  end
end

%% Wait for tasks to finish
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
  
  % Test files?
  torque_cleanup(ctrl);
end

end
