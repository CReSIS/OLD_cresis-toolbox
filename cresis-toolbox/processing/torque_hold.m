function torque_hold(ctrl,hold_status)
% torque_hold(ctrl,hold_status)
%
% Places holds on or removes holds from all the jobs specified.
%
% Inputs:
% ctrl = Several options which specify which batches to act on
%   1. Pass in a torque batch ctrl structure (only needs "batch_dir"
%      defined)
%   2. String containing regular expression which will be used to match
%      against the user for each job ('.*' cleans up all batches)
%   3. A vector of batch ids to apply hold to
% hold_status = mode must be one of the following
%   0: removes hold
%   1: applies hold
%   2: keyboard file
%   3: removes keyboard file
%
% Author: John Paden
%
% See also: torque_batch_list torque_cleanup torque_compile ...
%   torque_create_task torque_hold torque_job_list torque_job_status ...
%   torque_new_batch torque_print torque_rerun

%% Handle case where multiple batches may have been specified
if ~isstruct(ctrl)
  ctrls = torque_batch_list;
  for batch_idx = 1:length(ctrls)
    if (ischar(ctrl) && ~isempty(regexpi(ctrls{batch_idx}.user,ctrl))) ...
        || any(ctrls{batch_idx}.batch_id == ctrl)
      if hold_status == 2
        fprintf(' Placing keyboard on batch %d\n', ctrls{batch_idx}.batch_id);
      elseif hold_status == 3
        fprintf(' Removing keyboard on batch %d\n', ctrls{batch_idx}.batch_id);
      elseif hold_status == 1
        fprintf(' Placing user hold on jobs in batch %d\n', ctrls{batch_idx}.batch_id);
      elseif hold_status == 0
        fprintf(' Removing user hold on jobs in batch %d\n', ctrls{batch_idx}.batch_id);
      end
      torque_hold(ctrls{batch_idx},hold_status);
    else
      fprintf(' Skipping %d\n', ctrls{batch_idx}.batch_id);
    end
  end
  return;
end

%% This section actually does the placement/removal of holds

if hold_status == 2
  cmd = sprintf('touch %s', fullfile(ctrl.batch_dir,'keyboard'));
  try
    [status,result] = system(cmd);
  catch
    cmd
    warning('system call failed');
    %   keyboard;
  end
  return
  
elseif hold_status == 3
  keyboard_fn = fullfile(ctrl.batch_dir,'keyboard');
  if exist(keyboard_fn,'file')
    delete(keyboard_fn);
  end
  return
end

ctrl = torque_job_list(ctrl,ctrl.batch_id);

% For each job in the batch, remove/place as specified
for job_id = 1:length(ctrl.job_id_list)    
  if hold_status == 1
    if ctrl.job_status(job_id) == 'Q'
      cmd = sprintf('qhold -h u %i', ctrl.job_id_list(job_id));
      try
        [status,result] = system(cmd);
      catch
        cmd
        warning('system call failed');
        %   keyboard;
      end
    end
    
  elseif hold_status == 0
    if ctrl.job_status(job_id) == 'H'
      cmd = sprintf('qalter -h n %i', ctrl.job_id_list(job_id));
      try
        [status,result] = system(cmd);
      catch
        cmd
        warning('system call failed');
        %   keyboard;
      end
    end
  end
end

end
