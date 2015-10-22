function out = torque_submit_batch(fun,argsin)
% out = torque_submit_batch(fun,argsin)
%
% Runs an arbitrary job called "fun" with inputs "argsin". Introduces
% about 30 seconds of overhead as long as a compile is not necessary.
%
% INPUTS:
% fun = string containing the function name (e.g. 'hanning')
% argsin = cell vector of input arguments (e.g. {10}
% OUTPUTS:
% out = cell vector of output arguments
%
% out = torque_submit_batch('hanning',{10})
%
% Authors: John Paden
%
% See also: torque_batch_list torque_cleanup torque_compile ...
%   torque_create_task torque_hold torque_job_list torque_job_status ...
%   torque_new_batch torque_print torque_rerun

torque_compile(fun,[],0);

ctrl = torque_new_batch;
fprintf('Submitting %s\n', ctrl.batch_dir);

ctrl = torque_create_task(ctrl,str2func(fun),1,argsin);

ctrl = torque_rerun(ctrl);

[in,out] = torque_print(ctrl.batch_id,1,0);

torque_cleanup(ctrl);

return;


