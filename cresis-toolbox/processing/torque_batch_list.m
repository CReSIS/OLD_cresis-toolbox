function ctrls = torque_batch_list(param)
% ctrls = torque_batch_list(param)
%
% Gets a list of batches using the specified data location.
%
% Inputs:
% param = optional input (just needs to contain param.sched.data_location)
%   default value is to use global gRadar
%
% Outputs:
% ctrls = cell vector of batch jobs (can be used with the other torque_*
%   functions). Job information is not retrieved, just the basics.
%   Use torque_job_list to get job information.
%  .sched = scheduler structure
%   .data_location
%  .batch_dir = path to directory containing batch information
%  .job_id_file = path to file containing torque job ids
%  .in_path_dir = input arguments directory
%  .out_path_dir = output arguments directory
%  .stdout_path_dir = stdout directory
%  .error_path_dir = error directory
%  .user_file = path of file containing the username who owns batch
%  .user = string containing username who owns batch
%
% Author: John Paden
%
% See also: torque_batch_list torque_cleanup torque_compile ...
%   torque_create_task torque_hold torque_job_list torque_job_status ...
%   torque_new_batch torque_print torque_rerun

if ~exist('param','var')
  global gRadar;
  param = gRadar;
end

batch_dirs = get_filenames(param.sched.data_location,'batch_','','',struct('type','d'));

ctrls = [];
if nargout == 0
  fprintf('Batch Username    Directory\n');
end
for batch_idx = 1:length(batch_dirs)
  ctrls{batch_idx}.sched = param.sched;
  ctrls{batch_idx}.batch_dir = batch_dirs{batch_idx};
  ctrls{batch_idx}.job_id_file = fullfile(ctrls{batch_idx}.batch_dir,'job_id_file');

  [tmp batch_dir_name] = fileparts(ctrls{batch_idx}.batch_dir);
  ctrls{batch_idx}.batch_id = strtok(batch_dir_name(7:end),'_');
  ctrls{batch_idx}.batch_id = str2double(ctrls{batch_idx}.batch_id);
  ctrls{batch_idx}.in_path_dir = fullfile(ctrls{batch_idx}.batch_dir,'in');
  ctrls{batch_idx}.out_path_dir = fullfile(ctrls{batch_idx}.batch_dir,'out');
  ctrls{batch_idx}.stdout_path_dir = fullfile(ctrls{batch_idx}.batch_dir,'stdout');
  ctrls{batch_idx}.error_path_dir = fullfile(ctrls{batch_idx}.batch_dir,'error');

  ctrls{batch_idx}.user_file = fullfile(ctrls{batch_idx}.batch_dir,'user_file');
  [fid,msg] = fopen(ctrls{batch_idx}.user_file,'r');
  if fid < 1
    warning ('Could not open user file %s\n', ctrls{batch_idx}.user_file);
    ctrls{batch_idx}.user = '';
  else
    ctrls{batch_idx}.user = textscan(fid,'%s');
    fclose(fid);
    if isempty(ctrls{batch_idx}.user{1})
      warning ('User file is empty %s\n', ctrls{batch_idx}.user_file);
      ctrls{batch_idx}.user = '';
    else
      ctrls{batch_idx}.user = ctrls{batch_idx}.user{1}{1};
    end
  end
  
  if nargout == 0
    fprintf('%6d %12s %s\n', ctrls{batch_idx}.batch_id, ctrls{batch_idx}.user, ctrls{batch_idx}.batch_dir);
  end
end

end
