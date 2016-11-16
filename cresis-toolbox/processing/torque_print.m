function [in,out] = torque_print(batch,job_id,print_flag)
% [in,out] = torque_print(batch,job_id,print_flag)
%
% Prints out status information for a particular job ID.
%
% Inputs:
% batch: integer for which batch to get jobs info for
%   OR
%   control structure for the batch (e.g. returned from torque_job_list)
% job_id: OPTIONS:
%   integer for which job to print information about
%   string containing torque ID to print information about
% print_flag = optional (default is 1)
%  0: do not print status (just return output arguments)
%  1: print detailed status
%  2: print brief status
%
% Author: John Paden
%
% See also: torque_batch_list torque_cleanup torque_compile ...
%   torque_create_task torque_hold torque_job_list torque_job_status ...
%   torque_new_batch torque_print torque_rerun

if ~exist('print_flag','var') || isempty(print_flag)
  print_flag = 1;
end

if isstruct(batch)
  batch_id = batch.batch_id;
else
  batch_id = batch;
end

global gRadar;
ctrl.sched = gRadar.sched;

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
  error('Batch %d not found\n', batch_id);
end

[fid,msg] = fopen(ctrl.job_id_file,'r');
if fid < 1
  error('Could not open job id list file %s\n', ctrl.job_id_file);
  ctrl.job_id_list = [];
  ctrl.error_mask = [];
  ctrl.job_status = [];
  return;
end
ctrl.job_id_list = textscan(fid,'%f');
fclose(fid);
ctrl.job_id_list = ctrl.job_id_list{1};

if (print_flag ~=2) && (numel(job_id)>1)
  error('job_id must be a scalar for print_flag=0 and print_flag 0=1');
end

if print_flag == 0
  % Check for existence of input file
  fn = fullfile(ctrl.in_path_dir,sprintf('in_%d.mat',job_id));
  if exist(fn,'file')
    in = load(fn);
  else
    warning('No input file', fn);
  end
  
  % Check for existence of output file
  fn = fullfile(ctrl.out_path_dir,sprintf('out_%d.mat',job_id));
  out = load(fn);
  if exist(fn,'file')
    out = load(fn);
    if isfield(out,'errorstruct')
      fprintf('%s: %s\n', out.errorstruct.identifier, out.errorstruct.message);
      for stack_idx = 1:length(out.errorstruct.stack)
        fprintf('  %s: %d\n', out.errorstruct.stack(stack_idx).name, out.errorstruct.stack(stack_idx).line);
      end
    end
  else
    warning('No output file', fn);
  end
  return;
end

%%
if print_flag == 1
  if isnumeric(job_id)
    % Matlab side job ID
    torque_id = ctrl.job_id_list(job_id);
  else
    % Torque side job ID
    torque_id = str2double(job_id);
    job_id = ctrl.job_id_list(ctrl.job_id_list == torque_id)
  end
  
  
  cmd = sprintf('qstat -f %d', torque_id);
  status = -1;
  torque_attempts = 0;
  while status ~= 0
    try
      [status,result] = system(cmd);
    catch
      cmd
      warning('system call failedresources');
      keyboard;
    end
    if status ~= 0
      
      warning('qstat failed %d %s', status, result);
      torque_attempts = torque_attempts + 1;
      if torque_attempts >= 1
        break;
      end
      pause(2);
    end
  end
  
  % Print job status
  fprintf('Matlab Job ID %d\n', job_id);
  fprintf('Torque Job ID %d\n', torque_id);
  
  if status == 0
    fprintf('\n\n%s ======================================================\n', cmd);
    result
  end
  
  % Print stdout file
  fprintf('\n\nSTDOUT ======================================================\n');
  actual_job_id = find(ctrl.job_id_list==ctrl.job_id_list(job_id),1,'last');
  fn = fullfile(ctrl.stdout_path_dir,sprintf('stdout_%d.txt',actual_job_id));
  fprintf('%s\n', fn);
  if exist(fn,'file')
    type(fn);
  else
    fprintf('  Does not exist\n');
  end
  fprintf('\n\nERROR ======================================================\n');
  fn = fullfile(ctrl.error_path_dir,sprintf('error_%d.txt',actual_job_id));
  fprintf('%s\n', fn);
  if exist(fn,'file')
    type(fn);
  else
    fprintf('  Does not exist\n');
  end
  
  % Check for existence of input file
  fprintf('\nINPUT ======================================================\n');
  fn = fullfile(ctrl.in_path_dir,sprintf('in_%d.mat',job_id));
  fprintf('%s\n', fn);
  if exist(fn,'file')
    fprintf('  Exists\n');
    in = load(fn);
  else
    fprintf('  Does not exist\n');
    in = [];
  end
  
  % Check for existence of output file
  fprintf('\nOUTPUT ======================================================\n');
  fn = fullfile(ctrl.out_path_dir,sprintf('out_%d.mat',job_id));
  fprintf('%s\n', fn);
  if exist(fn,'file')
    fprintf('  Exists\n');
    out = load(fn);
    if isfield(out,'errorstruct')
      fprintf('%s: %s\n', out.errorstruct.identifier, out.errorstruct.message);
      for stack_idx = 1:length(out.errorstruct.stack)
        fprintf('  %s: %d\n', out.errorstruct.stack(stack_idx).name, out.errorstruct.stack(stack_idx).line);
      end
    end
    % Print job status
    fprintf('Matlab Job ID %d\n', job_id);
    fprintf('Torque Job ID %d\n', torque_id);
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
  
  if isempty(job_id)
    job_ids =ctrl.job_id_list;
    job_id = find(job_ids);
  end
  
  if isnumeric(job_id)
    % Matlab side job ID
    torque_id = ctrl.job_id_list(job_id);
  else
    % Torque side job ID
    torque_id = str2double(job_id);
    job_id = ctrl.job_id_list(ctrl.job_id_list == torque_id);
  end
  
  cmd = sprintf('qstat -f');
  
  status = -1;
  while status ~= 0
    try
      [status,result] = system(cmd);
    catch
      cmd
      warning('system call failedresources');
      keyboard;
    end
  end
  
  if status == 0
    qstat_res = textscan(result,'%s %s %s %s %s %s','HeaderLines',2,'Delimiter',sprintf(' \t'),'MultipleDelimsAsOne',1);
    [property, value] = qstat_res{1:2:3};
  end
  
  host_property_idx = find(strcmp(property,'exec_host'));
  mem_property_idx = find(strcmp(property,'resources_used.mem'));
  used_walltime_property_idx = find(strcmp(property,'resources_used.walltime'));
  pmem_property_idx = find(strcmp(property,'Resource_List.pmem'));
  list_walltime_property_idx = find(strcmp(property,'Resource_List.walltime'));
  
  num_jobes = numel(job_id);
  for i = 1 : num_jobes
    % Print job status
    fprintf('\nMatlab Job ID %d\n', job_id(i));
    fprintf('Torque Job ID %d\n', torque_id(i));
    
    
    host_property_value = value(host_property_idx(i));
    mem_property_value = value(mem_property_idx(i));
    used_walltime_property_value = value(used_walltime_property_idx(i));
    pmem_property_value = value(pmem_property_idx(i));
    list_walltime_property_value = value(list_walltime_property_idx(i));
    
    fprintf('\n%s: %s\n','exec_host',host_property_value{1});
    fprintf('%s: %s\n','resources_used.mem',mem_property_value{1});
    fprintf('%s: %s\n','resources_used.walltime',used_walltime_property_value{1});
    fprintf('%s: %s\n','Resource_List.pmem',pmem_property_value{1});
    fprintf('%s: %s\n','Resource_List.walltime',list_walltime_property_value{1});
    
    % Check for existence of input file
    fn = fullfile(ctrl.in_path_dir,sprintf('in_%d.mat',job_id(i)));
    if exist(fn,'file')
      fprintf('\nINPUT: Exists\n');
    else
      fprintf('\nINPUT: Does not exist\n');
    end
    
    % Check for existence of output file
    fn = fullfile(ctrl.out_path_dir,sprintf('out_%d.mat',job_id(i)));
    if exist(fn,'file')
      fprintf('OUTPUT: Exists\n');
    else
      fprintf('OUTPUT: Does not exist\n');
    end
    if i<num_jobes
      fprintf('==================================\n');
    end
  end
end
return

%% this is for a brief status of a single job 
% Print job statusjob_ids =numel(job_ids)
%   fprintf('Matlab Job ID %d\n', job_id);
%   fprintf('Torque Job ID %d\n', torque_id);
%   if print_flag == 2
%     if status == 0
%       qstat_res = textscan(result,'%s %s %s %s %s %s','HeaderLines',2,'Delimiter',sprintf(' \t'),'MultipleDelimsAsOne',1);
%       [property, value] = qstat_res{1:2:3};
%
%       host_property_idx = find(strcmp(property,'exec_host'),1,'last');
%       mem_property_idx = find(strcmp(property,'resources_used.mem'),1,'last');
%       used_walltime_property_idx = find(strcmp(property,'resources_used.walltime'),1,'last');
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
%       fn = fullfile(ctrl.in_path_dir,sprintf('in_%d.mat',job_id));
%       if exist(fn,'file')
%         fprintf('\nINPUT: Exists\n');
%       else
%         fprintf('\nINPUT: Does not exist\n');
%       end
%
%    % Check for existence of output file
%       fn = fullfile(ctrl.out_path_dir,sprintf('out_%d.mat',job_id));
%       if exist(fn,'file')
%         fprintf('OUTPUT: Exists\n');
%       else
%         fprintf('OUTPUT: Does not exist\n');
%       end
%     end
%   end
