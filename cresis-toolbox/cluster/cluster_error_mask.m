function cluster_error_mask(error_mask)
% cluster_error_mask(error_mask)
%
% Either creates all the variables needed to interpret the cluster error
% mask field OR it interprets the input error mask
%
% Inputs:
% error_mask: optional scalar integer
%
% Outputs:
% if no inputs, then the calling workspace will have variables added to it
%
% Author: John Paden
%
% See also: cluster_chain_stage, cluster_cleanup, cluster_compile
%   cluster_exec_job, cluster_get_batch, cluster_get_batch_list, 
%   cluster_hold, cluster_job, cluster_new_batch, cluster_new_task,
%   cluster_print, cluster_run, cluster_submit_batch, cluster_submit_task,
%   cluster_update_batch, cluster_update_task

if nargin == 0
  assignin('caller','out_fn_exist_error',1);
  assignin('caller','out_fn_load_error',2);
  assignin('caller','argsout_exist_error',4);
  assignin('caller','argsout_length_error',8);
  assignin('caller','errorstruct_exist_error',16);
  assignin('caller','errorstruct_contains_error',32);
  assignin('caller','success_error',64);
  assignin('caller','cluster_killed_error',128);
  assignin('caller','walltime_exceeded_error',256);
  assignin('caller','success_eval_error',512);
  assignin('caller','file_success_error',1024);
  assignin('caller','file_success_corrupt_error',2048);
  assignin('caller','max_mem_exceeded_error',4096);
  assignin('caller','insufficient_mcr_space',8192);
  
elseif error_mask == 0
  fprintf('No errors\n');
  
else
  cluster_error_mask;
  
  if bitand(error_mask,out_fn_exist_error)
    fprintf('out_fn_exist_error: Task output status file does not exist.\n'); 
  end
  if bitand(error_mask,out_fn_load_error)
    fprintf('out_fn_load_error: Task output status file is failing to load.\n'); 
  end
  if bitand(error_mask,argsout_exist_error)
    fprintf('argsout_exist_error: Task output status file does not contain argsout variable.\n'); 
  end
  if bitand(error_mask,argsout_length_error)
    fprintf('argsout_length_error: Task output status file argsout variable is the wrong length.\n'); 
  end
  if bitand(error_mask,errorstruct_exist_error)
    fprintf('errorstruct_exist_error: Task output status file does not contain errorstruct variable\n'); 
  end
  if bitand(error_mask,errorstruct_contains_error)
    fprintf('errorstruct_contains_error: Task output status file''s errorstruct contains an error. Task threw uncaught exception.\n'); 
  end
  if bitand(error_mask,success_error)
    fprintf('success_error: Success eval string returned that there was an error.\n'); 
  end
  if bitand(error_mask,cluster_killed_error)
    fprintf('cluster_killed_error: Cluster killed the job while running this task.\n'); 
  end
  if bitand(error_mask,walltime_exceeded_error)
    fprintf('walltime_exceeded_error: Cluster job exceeded the wall time.\n'); 
  end
  if bitand(error_mask,success_eval_error)
    fprintf('success_eval_error: Success eval string failed to evaluated and threw error.\n'); 
  end
  if bitand(error_mask,file_success_error)
    fprintf('file_success_error: One of the required output files was not created.\n'); 
  end
  if bitand(error_mask,file_success_corrupt_error)
    fprintf('file_success_corrupt_error: One of the required output files is corrupt.\n'); 
  end
  if bitand(error_mask,max_mem_exceeded_error)
    fprintf('max_mem_exceeded_error: The maximum memory was exceeded while running this job.\n'); 
  end
  if bitand(error_mask,insufficient_mcr_space)
    fprintf('insufficient_mcr_space: Job failed because there was insufficient space for the Matlab Compiler Runtime files.\n'); 
  end
end
