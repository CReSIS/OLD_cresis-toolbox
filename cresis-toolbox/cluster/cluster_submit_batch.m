function out = cluster_submit_batch(fun,argsin,num_args_out,cpu_time)
% out = cluster_submit_batch(fun,argsin,num_args_out,notes,cpu_time)
%
% Runs an arbitrary job called "fun" with inputs "argsin". Introduces
% about 30 seconds of overhead as long as a compile is not necessary.
%
% INPUTS:
% fun: string containing the function name (e.g. 'hanning') that the task
%   will run
% argsin: cell vector of input arguments to the task (e.g. {10,'arg2'})
% num_args_out: number of output arguments
% notes: string containing notes about the task
% cpu_time: expected maximum cpu time in seconds
% mem: expected maximum memory usage in bytes
%
% OUTPUTS:
% out: cell vector of output arguments from the task
%
% EXAMPLES:
% out = cluster_submit_batch('hanning',{10},1,60)
%
% Authors: John Paden
%
% See also: cluster_batch_list cluster_cleanup cluster_compile ...
%   cluster_create_task cluster_hold cluster_job_list cluster_job_status ...
%   cluster_new_batch cluster_print cluster_rerun

ctrl = cluster_new_batch;

cluster_compile(fun,[],0,ctrl);

param.task_function = fun;
param.argsin = argsin;
param.num_args_out = num_args_out;
param.cpu_time = cpu_time;
ctrl = cluster_create_task(ctrl,param,[]);

fprintf('Submitting %s\n', ctrl.batch_dir);

ctrl_chain = {{ctrl}};

ctrl_chain = cluster_run(ctrl_chain);

[in,out] = cluster_print(ctrl_chain{1}{1}.batch_id,1,0);

cluster_cleanup(ctrl);

return;
