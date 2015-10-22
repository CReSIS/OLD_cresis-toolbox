% script demo_parallel_proc
%
% This shows the basics of using the parallel processing
% toolbox on a local machine with multiple CPUS/cores.
%
% drebber has 8 processors
% small has X processors

%% Cluster setup

% param.sched_type = scheduler type
%   (local sets up a local scheduler)
param.sched_type = 'local';

% local_sched = local scheduler
local_sched = findResource('scheduler','type',param.sched_type);

%% Example job creation, adding tasks, submitting/execution,
% waiting for completion, getting arguments, and clean up

% Construct a job object using the default configuration.
job = createJob(local_sched);

% Add 10 tasks to the job.
%  - Tasks are defined by a function handle
%  - The function handle here is "@rand" --> Matlab's rand function
%  - Function handles are usually identified by the "@" or str2func
%    commands
for task_ind = 1:10
  createTask(job, @rand, 1, {10});
end

% Run the job.
fprintf('Submitting job %s %d\n', get(job,'Name'), get(job,'ID'));
submit(job);

% Wait for the job to finish running, and retrieve the job results.
waitForState(job);
out = getAllOutputArguments(job);

% Display the random matrix returned from the third task.
disp(out{3});

% Destroy the job
destroy(job);

return;
