function create_task_debug(ctrl)
% create_task_debug(ctrl)
%
% This function can be run when the "ctrl" variable exists that create_task.m
% uses. It prints out debug information for each job
%
% Recommendation:
%  global gCtrl;
%  gCtrl = ctrl;
%  create_task_debug(gCtrl);
%
% Author: John Paden

for job_idx = 1:length(ctrl.jobs)
  if ishandle(ctrl.jobs{job_idx}.job)
    get(ctrl.jobs{job_idx}.job)
  end
end

return;
