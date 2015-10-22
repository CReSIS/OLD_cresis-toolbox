function torque_cleanup(ctrl)
% torque_cleanup(ctrl)
%
% Delete batches (jobs from queue and temporary files)
%
% Inputs:
% ctrl = Several options which specify which batches to act on
%   1. Pass in a torque batch ctrl structure (only needs "batch_dir"
%      defined)
%   2. String containing regular expression which will be used to match
%      against the user for each job ('.*' cleans up all batches)
%   3. A vector of batch ids to apply hold to
%
% Examples:
%   torque_cleanup(ctrl)
%   torque_cleanup('user_name')
%   torque_cleanup([1 3])
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
      fprintf(' Deleting jobs in batch %d\n', ctrls{batch_idx}.batch_id);
      torque_cleanup(ctrls{batch_idx});
    else
      fprintf(' Skipping %d\n', ctrls{batch_idx}.batch_id);
    end
  end
  return;
end

%% This section actually does the placement/removal of holds
ctrl = torque_job_list(ctrl,ctrl.batch_id);

% For each job in the batch, delete the job
for job_id = 1:length(ctrl.job_id_list)
  if ctrl.job_status(job_id) ~= 'C'
    % Only delete jobs that have not been completed (completed jobs
    % are effectively deleted already)
    cmd = sprintf('qdel %i', ctrl.job_id_list(job_id));
    try
      [status,result] = system(cmd);
    catch
      cmd
      warning('system call failed');
      %   keyboard;
    end
  end
end

% Finally, remove the batch directory containing all the jobs' information
if exist(ctrl.batch_dir,'dir')
  fprintf('  Removing %s\n', ctrl.batch_dir);
  
  status = -1;
  torque_attempts = 0;
  while status ~= 0
    try
      if exist(ctrl.batch_dir,'dir')
        rmdir(ctrl.batch_dir,'s');
      end
      status = 0;
    catch ME
      warning('rmdir failed');
      torque_attempts = torque_attempts + 1;
      pause(3);
      if torque_attempts > 2
        % There is potentially something wrong with the torque tasks if
        % rmdir fails this many times. Look into it before running
        % "dbcont" to keep trying to remove the directory.
        keyboard;
      end
    end
  end
end

end
