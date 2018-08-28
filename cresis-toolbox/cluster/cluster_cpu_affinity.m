function cluster_cpu_affinity(num_my_processors,process_list)
% cluster_cpu_affinity(num_my_processors,process_list)
%
% Assigns the current process to "num_my_processors" using cpu affinity
% settings.
%
% Inputs:
% num_my_processors: the number of physical cores to assign.
% process_list: string argument which overrides process list to taskset
%
% Author: John Paden
%
% See also: cluster_chain_stage, cluster_cleanup, cluster_compile
%   cluster_exec_job, cluster_get_batch, cluster_get_batch_list,
%   cluster_hold, cluster_job, cluster_new_batch, cluster_new_task,
%   cluster_print, cluster_run, cluster_submit_batch, cluster_submit_task,
%   cluster_update_batch, cluster_update_task

if ~exist('process_list','var')
  process_list = '';
end

try
  % 1. Get Matlab's process id:
  % feature getpid
  pid = feature('getpid');
  
  % 2. Determine what processors this machine has
  % numactl -H
  cmd = 'numactl -H';
  [status,result] = system(cmd);
  node_idxs = regexpi(result,'node [0-9] cpus:');
  node = [];
  for node_idxs_idx = 1:length(node_idxs)
    end_idx = find(result(node_idxs(node_idxs_idx)+12:end) == 10,1);
    node(end+1,:) = cellfun(@str2double, regexpi(result(node_idxs(node_idxs_idx)+12+(1:end_idx-2)),' ','split'));
  end
  
  % 3. Determine what the cpu affinity is for each process
  cmd = 'for pid in $(ps -a -o pid=); do taskset -p $pid 2>/dev/null | awk ''{print $6}''; done';
  [status,result] = system(cmd);
  
  process = [];
  if status == 0
    cpu_affinity_dec = sscanf(result, '%x');
    for proc_idx = 1:length(cpu_affinity_dec)
      new_process = double(dec2bin(cpu_affinity_dec(proc_idx),numel(node)) - '0').';
      if ~all(new_process)
        process(end+1,:) = new_process;
      end
    end
  else
    error('Failed to get cpu affinity.');
  end
  
  % Identify if any processes have affinity. Avoid processors that have
  % affinity to specific processes.
  if isempty(process)
    % No processes have affinity, so assume that all processors are available
    processors = zeros(1,numel(node));
  else
    processors = fliplr(any(process,1));
  end
  
  % Loop through available processors and select processors that are
  % available and are not logical pairs to already selected processors. If
  % there are not enough processors available, print a warning message since
  % this never happen.Confirm which process is the lead logical pr
  my_processors = [];
  for cpu_idx = 1:length(processors)
    if processors(cpu_idx) == 0
      cpu = cpu_idx - 1;
      % Find all of this processor's logical pairs
      cmd = sprintf('cat /sys/devices/system/cpu/cpu%d/topology/thread_siblings_list', cpu);
      [status,result] = system(cmd);
      logical_processors = cellfun(@str2double,regexp(result,',','split'));
      % Verify that none of the processor's logical pairs are used
      if ~all(~processors(logical_processors+1))
        processors(logical_processors+1) = 0;
        continue;
      end
      % Take this processor and its logical pairs
      processors(logical_processors+1) = 1;
      my_processors(end+1) = cpu;
      num_my_processors = num_my_processors - 1;
    end
    if num_my_processors == 0
      break;
    end
  end
  if num_my_processors ~= 0
    error('Insufficient processors available.');
  end
  
  % 4. Set this processes CPU affinity
  % taskset -p MY_CPU_AFFINITY PID
  process = zeros(1,size(process,2));
  process(my_processors+1) = 1;
  process_str = sprintf('%d,',find(process)-1);
  process_str = process_str(1:end-1);
  
  % Set parent process CPU affinity
  cmd = sprintf('ps -p %d -o pid,ppid | grep %d | awk ''{print $2}''', pid, pid);
  [status,result] = system(cmd);
  ppid = cellfun(@str2double,regexp(result,'\n','split'));
  if ~isempty(process_list)
    % Override
    cmd = sprintf('taskset -cp %s %d', process_list, ppid(1));
  else
    cmd = sprintf('taskset -cp %s %d', process_str, ppid(1));
  end
  fprintf('%s\n', cmd);
  system(cmd);
  
  % Set this process and all of its thread's cpu affinity
  cmd = sprintf('ps -eLo pid,tid | grep %d| awk ''{print $2}''', pid);
  [status,result] = system(cmd);
  pids = cellfun(@str2double,regexp(result,'\n','split'));
  for pid_idx = 1:length(pids)
    pid = pids(pid_idx);
    if ~isnan(pid)
      if ~isempty(process_list)
        % Override
        cmd = sprintf('taskset -cp %s %d >/dev/null 2>&1', process_list, pid);
      else
        cmd = sprintf('taskset -cp %s %d >/dev/null 2>&1', process_str, pid);
      end
      % fprintf('%s\n', cmd);
      system(cmd);
    end
  end
  
catch ME
  warning(ME.getReport);
end
