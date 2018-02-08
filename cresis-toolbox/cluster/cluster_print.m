function [in,out] = cluster_print(ctrl,task_id,print_flag)
% [in,out] = cluster_print(ctrl,task_id,print_flag)
%
% Prints out status information for a particular job ID.
%
% Inputs:
% ctrl: integer for which batch to get jobs info for
%   OR
%   control structure for the batch (e.g. returned from cluster_job_list)
% task_id: OPTIONS:
%   integer for which job to print information about
%   string containing torque ID to print information about
% print_flag = optional (default is 1)
%  0: do not print status (just return output arguments)
%  1: print detailed status
%  2: print brief status
%
% Author: John Paden
%
% See also: cluster_chain_stage, cluster_cleanup, cluster_compile
%   cluster_exec_job, cluster_get_batch, cluster_get_batch_list, 
%   cluster_hold, cluster_job, cluster_new_batch, cluster_new_task,
%   cluster_print, cluster_run, cluster_submit_batch, cluster_submit_task,
%   cluster_update_batch, cluster_update_task

if ~exist('print_flag','var') || isempty(print_flag)
  print_flag = 1;
end

if isstruct(ctrl)
  batch_id = ctrl.batch_id;
  ctrls = cluster_get_batch_list(ctrl);
else
  batch_id = ctrl;
  ctrls = cluster_get_batch_list;
end

found = 0;
global gRadar;
for batch_idx = 1:length(ctrls)
  if ctrls{batch_idx}.batch_id == batch_id
    found = 1;
    ctrl = ctrls{batch_idx};
    break;
  end
end

if found == 0
  error('Batch %d not found\n', batch_id);
end

[fid,msg] = fopen(ctrl.job_id_fn,'r');
if fid < 1
  error('Could not open job id list file %s\n', ctrl.job_id_fn);
  ctrl.job_id_list = [];
  
  ctrl.error_mask = [];
  ctrl.job_status = [];
  return;
end
ctrl.job_id_list = textscan(fid,'%f');
fclose(fid);
ctrl.job_id_list = ctrl.job_id_list{1};

if (print_flag ~=2) && (numel(task_id)>1)
  error('task_id must be a scalar for print_flag=0 and print_flag 0=1');
end

if print_flag == 0
  % Check for existence of input files
  
  % Create in/out filenames
  static_in_fn = fullfile(ctrl.in_fn_dir,'static.mat');
  dynamic_in_fn = fullfile(ctrl.in_fn_dir,'dynamic.mat');
  
  if exist(static_in_fn,'file') && exist(dynamic_in_fn,'file')
    % Read input
    sparam = load(static_in_fn);
    dparam_task_field = sprintf('dparam_%d',task_id);
    dparam = load(dynamic_in_fn,dparam_task_field);
    if ~isfield(dparam,dparam_task_field)
      error('Task input does not exist for task %d', task_id);
    end
    in = merge_structs(sparam.static_param,dparam.(dparam_task_field));
  else
    warning('Missing one of the input files:\n  %s\n  %s', static_in_fn, dynamic_in_fn);
    in = [];
  end
  
  % Check for existence of output file
  out_fn = fullfile(ctrl.out_fn_dir,sprintf('out_%d.mat',task_id));
  if exist(out_fn,'file')
    out = load(out_fn);
    if isfield(out,'errorstruct')
      if ~isempty(out.errorstruct)
        fprintf('%s: %s\n', out.errorstruct.identifier, out.errorstruct.message);
        for stack_idx = 1:length(out.errorstruct.stack)
          fprintf('  %s: %d\n', out.errorstruct.stack(stack_idx).name, out.errorstruct.stack(stack_idx).line);
        end
      end
    end
  else
    warning('No output file', out_fn);
    out = [];
  end
  return;
end

%%
if print_flag == 1
  if isnumeric(task_id)
    % Matlab side job ID
    job_id = ctrl.job_id_list(task_id);
  else
    % Torque side job ID
    job_id = str2double(task_id);
    task_id = ctrl.job_id_list(ctrl.job_id_list == job_id)
  end
  
  
  if strcmpi(ctrl.cluster.type,'torque')
    cmd = sprintf('qstat -f %d  </dev/null', job_id);
    try; [status,result] = system(cmd); end;
    
  elseif strcmpi(ctrl.cluster.type,'matlab')
    cmd = 'NA';
    status = -1;
    
  elseif strcmpi(ctrl.cluster.type,'slurm')
    %cmd = sprintf('sstat --format=AveCPU,AvePages,AveRSS,AveVMSize,JobID -j %d --allsteps', job_id);
    cmd = sprintf('scontrol show job %d', job_id);
    try; [status,result] = system(cmd); end;
    
  elseif strcmpi(ctrl.cluster.type,'debug')
    cmd = 'NA';
    status = -1;
  end
  
  % Print job status
  fprintf('Matlab Task ID %d\n', task_id);
  fprintf('Cluster Job ID %d\n', job_id);
  
  if status == 0
    fprintf('\n\n%s ======================================================\n', cmd);
    result
  end
  
  % Print stdout file
  if any(strcmpi(ctrl.cluster.type,{'torque','slurm'}))
    fprintf('\n\nSTDOUT ======================================================\n');
    actual_job_id = find(ctrl.job_id_list==ctrl.job_id_list(task_id),1,'last');
    fn = fullfile(ctrl.stdout_fn_dir,sprintf('stdout_%d.txt',actual_job_id));
    fprintf('%s\n', fn);
    if exist(fn,'file')
      type(fn);
    else
      fprintf('  Does not exist\n');
    end
    fprintf('\n\nERROR ======================================================\n');
    fn = fullfile(ctrl.error_fn_dir,sprintf('error_%d.txt',actual_job_id));
    fprintf('%s\n', fn);
    if exist(fn,'file')
      type(fn);
    else
      fprintf('  Does not exist\n');
    end
  elseif strcmpi(ctrl.cluster.type,'matlab')
    fprintf('\n\nSTDOUT ======================================================\n');
    ctrl.cluster.jm.Jobs.findobj('ID',ctrl.job_id_list(task_id)).Tasks.Diary
  end
  
  % Check for existence of input file
  fprintf('\nINPUT ======================================================\n');
  static_in_fn = fullfile(ctrl.in_fn_dir,'static.mat');
  fprintf('%s\n', static_in_fn);
  if exist(static_in_fn,'file')
    fprintf('  Exists\n');
    sparam = load(static_in_fn);
  else
    fprintf('  Does not exist\n');
    in = [];
  end
  dynamic_in_fn = fullfile(ctrl.in_fn_dir,'dynamic.mat');
  fprintf('%s\n', dynamic_in_fn);
  if exist(dynamic_in_fn,'file')
    fprintf('  Exists\n');
    dparam_task_field = sprintf('dparam_%d',task_id);
    dparam = load(dynamic_in_fn,dparam_task_field);
    in = merge_structs(sparam.static_param,dparam.(dparam_task_field));
  else
    fprintf('  Does not exist\n');
    in = [];
  end
  
  % Check for existence of output file
  fprintf('\nOUTPUT ======================================================\n');
  fn = fullfile(ctrl.out_fn_dir,sprintf('out_%d.mat',task_id));
  fprintf('%s\n', fn);
  if exist(fn,'file')
    fprintf('  Exists\n');
    out = load(fn);
    if isfield(out,'errorstruct') && ~isempty(out.errorstruct)
      fprintf('%s: %s\n', out.errorstruct.identifier, out.errorstruct.message);
      for stack_idx = 1:length(out.errorstruct.stack)
        fprintf('  %s: %d\n', out.errorstruct.stack(stack_idx).name, out.errorstruct.stack(stack_idx).line);
      end
    end
  else
    fprintf('  Does not exist\n');
    out = [];
  end
  
  if nargout == 0
    clear in out
  elseif nargout == 1
    clear out
  end
end

%% brief status of any set of jobes
if print_flag == 2
  if max(task_id)>numel(ctrl.job_id_list)
    Num_jobs = numel(ctrl.job_id_list);
    fprintf('\n%s: %d\n','Job ID not found',task_id(find(task_id>Num_jobs,1)))
    return;
  end
  if isnumeric(task_id)
    % Matlab side job ID
    job_id = ctrl.job_id_list(task_id);
  else
    % Torque side job ID
    job_id = str2double(task_id);
    task_id = ctrl.job_id_list(ctrl.job_id_list == job_id);
  end
  
  for idx = 1 : numel(job_id)
    
    fprintf('\nMatlab Job ID %d\n', task_id(idx));
    fprintf('Torque Job ID %d\n', job_id(idx));
    
    if strcmpi(ctrl.cluster.type,'custom_torque')
      cmd = sprintf('qstat -f %d  </dev/null', job_id);
      [status,result] = robust_system(cmd,0);
      
    elseif strcmpi(ctrl.cluster.type,'matlab')
      keyboard
      
    elseif strcmpi(ctrl.cluster.type,'debug')
      cmd = 'NA';
      status = 0;
      result = '';
    end
    
    if status==0 && isempty(result)
      qstat_res = textscan(result,'%s %s %s %s %s %s','HeaderLines',2,'Delimiter',sprintf(' \t'),'MultipleDelimsAsOne',1);
      [property, value] = qstat_res{[1 3]};
      
      % DEBUG: TORQUE ID
      host_property_idx = find(strcmp(property,'exec_host'));
      mem_property_idx = find(strcmp(property,'resources_used.mem'));
      used_walltime_property_idx = find(strcmp(property,'resources_used.walltime'));
      pmem_property_idx = find(strcmp(property,'Resource_List.pmem'));
      list_walltime_property_idx = find(strcmp(property,'Resource_List.walltime'));
      
      % Print job status
      
      % DEBUG: Find matching index in qstat -f result based on matching job_id(i)
      host_property_value = value(host_property_idx);
      mem_property_value = value(mem_property_idx);
      used_walltime_property_value = value(used_walltime_property_idx);
      pmem_property_value = value(pmem_property_idx);
      list_walltime_property_value = value(list_walltime_property_idx);
      
      if isempty(host_property_value)
        fprintf('\n%s: %s\n','exec_host','Not found');
      else
        fprintf('\n%s: %s\n','exec_host',host_property_value{1});
      end
      if isempty(mem_property_value)
        fprintf('%s: %s\n','resources_used.mem','Not found');
      else
        fprintf('%s: %s\n','resources_used.mem',mem_property_value{1});
      end
      if isempty(used_walltime_property_value)
        fprintf('%s: %s\n','resources_used.walltime','Not found');
      else
        fprintf('%s: %s\n','resources_used.walltime',used_walltime_property_value{1});
      end
      if isempty(pmem_property_value)
        fprintf('%s: %s\n','Resource_List.pmem','Not found');
      else
        fprintf('%s: %s\n','Resource_List.pmem',pmem_property_value{1});
      end
      if isempty(list_walltime_property_value)
        fprintf('%s: %s\n','Resource_List.walltime','Not found');
      else
        fprintf('%s: %s\n','Resource_List.walltime',list_walltime_property_value{1});
      end
      
      % Check for existence of input file
      fn = fullfile(ctrl.in_fn_dir,sprintf('in_%d.mat',task_id(idx)));
      if exist(fn,'file')
        fprintf('\nINPUT: Exists\n');
      else
        fprintf('\nINPUT: Does not exist\n');
      end
      
      % Check for existence of output file
      fn = fullfile(ctrl.out_fn_dir,sprintf('out_%d.mat',task_id(idx)));
      if exist(fn,'file')
        fprintf('OUTPUT: Exists\n');
      else
        fprintf('OUTPUT: Does not exist\n');
      end
      
    end
    fprintf('==================================\n');
  end
end
return

%% this is for a brief status of a single job
% Print job statusjob_ids =numel(job_ids)
%   fprintf('Matlab Job ID %d\n', task_id);
%   fprintf('Torque Job ID %d\n', job_id);
%   if print_flag == 2
%     if status == 0
%       qstat_res = textscan(result,'%s %s %s %s %s %s','HeaderLines',2,'Delimiter',sprintf(' \t'),'MultipleDelimsAsOne',1);
%       [property, value] = qstat_res{1:2:3};
%
%       host_property_idx = find(strcmp(property,'exec_host'),1,'last');
%       mem_property_idx = find(strcmp(property,'resources_used.mem'),1,'last');
%       used_walltime_property_idx = ,print_flag)
%find(strcmp(property,'resources_used.walltime'),1,'last');
%       pmem_property_idx = find(strcmp(property,'Resource_List.pmem'),1,'last');
%       list_walltime_property_idx = find(strcmp(property,'Resource_List.walltime'),1,'last');
%
%       host_property_value = value(host_property_idx);
%       mem_property_value = value(mem_property_idx);
%       used_walltime_property_value = value(used_walltime_property_idx);
%       pmem_property_value = value(pmem_property_idx);
%       list_walltime_property_value = value(list_walltime_property_idx);
%       if print_flag == 2
%       fprintf('\n%s: %s\n','exec_host',host_property_value{1});
%       fprintf('%s: %s\n','resources_used.mem',mem_property_value{1});
%       fprintf('%s: %s\n','resources_used.walltime',used_walltime_property_value{1});
%       fprintf('%s: %s\n','Resource_List.pmem',pmem_property_value{1});
%       fprintf('%s: %s\n','Resource_List.walltime',list_walltime_property_value{1});
%
%    % Check for existence of input file
%       fn = fullfile(ctrl.in_fn_dir,sprintf('in_%d.mat',task_id));
%       if exist(fn,'file')
%         fprintf('\nINPUT: Exists\n');
%       else
%         fprintf('\nINPUT: Does not exist\n');
%       end
%
%    % Check for existence of output file
%       fn = fullfile(ctrl.out_fn_dir,sprintf('out_%d.mat',task_id));
%       if exist(fn,'file')
%         fprintf('OUTPUT: Exists\n');
%       else
%         fprintf('OUTPUT: Does not exist\n');
%       end
%     end
%   end
