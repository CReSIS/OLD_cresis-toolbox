function ctrl = torque_job_list(ctrl,batch_id)
% ctrl = torque_job_list(ctrl,batch_id)
%
% Updates job status information from the torque cluster. Also prints
% out status information when job changes status.
%
% Inputs:
% ctrl = Default is to leave empty, []. However, if ctrl is provided,
%  the ctrl.sched.data_location variable is used to get the batch list
%  rather than gRadar.sched.data_location.
% batch_id: integer for which batch to get jobs info for
%
% Outputs:
% ctrl = ctrl structure for the specified batch_id
%
% C -  Job is completed after having run/
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

if isempty(ctrl)
  global gRadar;
  ctrl.sched = gRadar.sched;
end

ctrls = torque_batch_list(ctrl);

found = 0;
for batch_idx = 1:length(ctrls)
  if ctrls{batch_idx}.batch_id == batch_id
    found = 1;
    ctrl = ctrls{batch_idx};
    break;
  end
end

if found == 0
  fprintf('Batch %d not found\n', batch_id);
end

[fid,msg] = fopen(ctrl.job_id_file,'r');
if fid < 1
  warning ('Could not open job id list file %s\n', ctrl.job_id_file);
  ctrl.job_id_list = [];
  ctrl.error_mask = [];
  ctrl.job_status = [];
  return;
end
ctrl.job_id_list = textscan(fid,'%f');
fclose(fid);
ctrl.job_id_list = ctrl.job_id_list{1};
ctrl.job_id = length(ctrl.job_id_list);
ctrl.submission_queue = [];

cmd = 'qstat </dev/null';
try
  [status,result] = system(cmd);
catch
  cmd
  warning('system call failed');
  %   keyboard;
end

job_status_found = zeros(size(ctrl.job_id_list));
ctrl.error_mask = zeros(size(ctrl.job_id_list)).';
if status == 0
  if ~isempty(result)
    qstat_res = textscan(result,'%s %s %s %s %s %s','HeaderLines',2,'Delimiter',sprintf(' \t'),'MultipleDelimsAsOne',1);
    for idx = 1:size(qstat_res{1},1)
      qstat_res{7}(idx) = str2double(strtok(qstat_res{1}{idx},'.'));
    end
    for idx = 1:size(qstat_res{5},1)
      job_ids = find(qstat_res{7}(idx)==ctrl.job_id_list);
      % Qstat returns all jobs, just look at jobs in this batch
      while ~isempty(job_ids)
        job_id = job_ids(1);
        job_ids = job_ids(2:end);
        job_status_found(job_id) = 1;
        ctrl.job_status(job_id) = qstat_res{5}{idx};
        if ctrl.job_status(job_id) == 'C'
          in_path = fullfile(ctrl.in_path_dir,sprintf('in_%d.mat',job_id));
          clear param; load(in_path,'param');
          if param.conforming
            % Good output implies: output file exists and contains a first
            % output argument containing the number "1"
            out_path = fullfile(ctrl.out_path_dir,sprintf('out_%d.mat',job_id));
            if ~exist(out_path,'file')
              % Output file does not exist
              ctrl.error_mask(job_id) = 1;
            else
              try
                ctrl.out{job_id} = load(out_path,'argsout');
              catch
                ctrl.out{job_id} = [];
              end
              if isempty(ctrl.out{job_id})
                % Corrupt output file, missing argsout, or empty argsout
                ctrl.error_mask(job_id) = 2;
              else
                try
                  if ctrl.out{job_id}.argsout{1} ~= 1
                    % Failure in job's Matlab code (i.e. probably exitted
                    % due to an "error" exception in users code)
                    ctrl.error_mask(job_id) = 3;
                  end
                catch
                  % Non-conforming output (e.g. argsout{1} is not double scalar)
                  ctrl.error_mask(job_id) = 4;
                end
              end
            end
          else
            out_path = fullfile(ctrl.out_path_dir,sprintf('out_%d.mat',job_id));
            if ~exist(out_path,'file')
              % Output file does not exist
              ctrl.error_mask(job_id) = 1;
            else
              try
                ctrl.out{job_id} = load(out_path,'argsout');
              catch
                % Corrupt output file, missing argsout
                ctrl.error_mask(job_id) = 2;
                ctrl.out{job_id} = [];
              end
            end
          end
        end
      end
    end
  end
end

if any(job_status_found==0)
  lost_jobs = find(~job_status_found);
  for job_id = reshape(lost_jobs,[1 length(lost_jobs)])
    ctrl.job_status(job_id) = 'C';
    in_path = fullfile(ctrl.in_path_dir,sprintf('in_%d.mat',job_id));
    if exist(in_path,'file')
      clear param; load(in_path,'param');
    else
      warning('Input file %s missing. Assuming non-conforming.', in_path);
      clear param; param.conforming = false;
    end
    if param.conforming
      % Good output implies: output file exists and contains a first
      % output argument containing the number "1"
      out_path = fullfile(ctrl.out_path_dir,sprintf('out_%d.mat',job_id));
      if ~exist(out_path,'file')
        % Output file does not exist
        ctrl.error_mask(job_id) = 1;
      else
        try
          ctrl.out{job_id} = load(out_path,'argsout');
        catch
          ctrl.out{job_id} = [];
        end
        if isempty(ctrl.out{job_id})
          % Corrupt output file, missing argsout, or empty argsout
          ctrl.error_mask(job_id) = 2;
        else
          try
            if ctrl.out{job_id}.argsout{1} ~= 1
              % Failure in job's Matlab code (i.e. probably exitted
              % due to an "error" exception in users code)
              ctrl.error_mask(job_id) = 3;
            end
          catch
            % Non-conforming output (e.g. argsout{1} is not double scalar)
            ctrl.error_mask(job_id) = 4;
          end
        end
      end
    else
      out_path = fullfile(ctrl.out_path_dir,sprintf('out_%d.mat',job_id));
      if ~exist(out_path,'file')
        % Output file does not exist
        ctrl.error_mask(job_id) = 1;
      else
        try
          ctrl.out{job_id} = load(out_path,'argsout');
        catch
          % Corrupt output file, missing argsout
          ctrl.error_mask(job_id) = 2;
          ctrl.out{job_id} = [];
        end
      end
    end
  end
end

end
