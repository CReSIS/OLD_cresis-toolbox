% script demo_matlab_MDCS
%
% This shows the basics of using a Matlab cluster running
% Matlab's Distributed Computing Server.
%
% CReSIS has MDCS running on their cluster.
%   Currently 128 nodes split into two schedulers (jm and jm2)
%   with 64 nodes each.

%% Cluster setup

% param.sched_type = scheduler type
%   (jobmanager is Matlab default, torque is another option for later)
param.sched_type = 'jobmanager';
% param.sched_url = Matlab Distributed Computing Server (MDCS)
param.sched_url = 'heimdall.cluster.cresis.ku.edu';
% param.sched_name = schedule name (options are jm and jm2)
param.sched_name = 'jm';

% jm = job manager (jm and 1 of 2 running on heimdall)
jm = findResource('scheduler','type',param.sched_type, ...
  'LookupUrl',param.sched_url,'name',param.sched_name);

%% Example job creation, adding tasks, submitting/execution,
% waiting for completion, getting arguments, and clean up

% Construct a job object using the default configuration.
job = createJob(jm);

% Add 10 tasks to the job.
%  - Tasks are defined by a function handle
%  - The function handle here is "@rand" --> Matlab's rand function
%  - Function handles are usually identified by the "@" or str2func
%    commands
for task_ind = 1:10
  createTask(job, @rand, 1, {10},'MaximumNumberOfRetries',2);
end

% Run the job.
submit(job);

% Wait for the job to finish running, and retrieve the job results.
waitForState(job);
out = getAllOutputArguments(job);

% Display the random matrix returned from the third task.
disp(out{1});

% Destroy the job
destroy(job);

return;
