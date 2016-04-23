function ctrl = torque_check_conforming(ctrl,job_id)
% ctrl = torque_check_conforming(ctrl,job_id)
%
% Updates a specific job status. Support function for torque_job_status.
%
% Inputs:
% ctrl: ctrl structure returned from torque_new_batch
%  .job_id_list = Nx1 vector of torque job IDs
%  .job_status = Nx1 vector of job status
%  .out_path_dir = string containing the output directory
%  .error_mask = Nx1 vector of error status
%  .out = up-to-Nx1 vector of cells containing the outputs for each
%    job as they complete (read in from the output.mat files)
% job_id: the job id to check up on
%
% Outputs:
% ctrl = updated ctrl structure
%
% C -  Job is completed after having run.
% E -  Job is exiting after having run.
% H -  Job is held.
% Q -  job is queued, eligible to run or routed.
% R -  job is running.
% T -  job is being moved to new location.
% W -  job is waiting for its execution time
%      (-a option) to be reached.
% S -  (Unicos only) job is suspend.
%
% Author: John Paden
%
% See also: torque_batch_list torque_cleanup torque_compile ...
%   torque_create_task torque_hold torque_job_list torque_job_status ...
%   torque_new_batch torque_print torque_rerun

% Good output implies: output file exists and contains a first
% output argument containing the number "1"
out_path = fullfile(ctrl.out_path_dir,sprintf('out_%d.mat',job_id));

% If job state is 'C':
%  If output does not exist or is in error, then we report an error
% If job status is anything else (e.g. 'E' or 'R'):
%  If output does not exist or is in error, then we wait and try again
%  later

if ~exist(out_path,'file')
  % Output file does not exist
  if ctrl.job_status(job_id) == 'C'
    ctrl.error_mask(job_id) = 1;
  end
  
else
  try
    ctrl.out{job_id} = load(out_path,'argsout');
  catch
    ctrl.out{job_id} = [];
  end
  if isempty(ctrl.out{job_id})
    % Corrupt output file
    if ctrl.job_status(job_id) == 'C'
      ctrl.error_mask(job_id) = 2;
    end
  else
    try
      if ~isfield(ctrl.out{job_id},'argsout') || isempty(ctrl.out{job_id}.argsout)
        if ctrl.job_status(job_id) == 'C'
          % Failure in job's Matlab code (i.e. probably exitted
          % due to an "error" exception in users code)
          ctrl.error_mask(job_id) = 3;
        end
      elseif ctrl.out{job_id}.argsout{1} ~= 1
        % Success!
        ctrl.error_mask(job_id) = 0;
        ctrl.job_status(job_id) = 'C';
      end
    catch
      % Non-conforming output (e.g. argsout{1} is not double scalar)
      ctrl.error_mask(job_id) = 4;
      ctrl.job_status(job_id) = 'C';
    end
  end
  
  if ctrl.error_mask(job_id) && ctrl.job_status(job_id) == 'C'
    if exist(out_path,'file')
      out = load(out_path);
      error_string = '';
      if isfield(out,'errorstruct')
        error_string = sprintf('%s: %s\n', out.errorstruct.identifier, out.errorstruct.message);
        for stack_idx = 1:length(out.errorstruct.stack)
          error_string = cat(2,error_string,...
            sprintf('  %s: %d\n', out.errorstruct.stack(stack_idx).name, out.errorstruct.stack(stack_idx).line));
        end
      end
      warning('Job %d:%d/%d Error:\n%s', ctrl.batch_id, job_id, ctrl.job_id_list(job_id), error_string);
    else
      warning('Job %d:%d/%d Error, but no output file error message\n', ctrl.batch_id, job_id, ctrl.job_id_list(job_id));
      out = [];
    end
  end
end

return;
