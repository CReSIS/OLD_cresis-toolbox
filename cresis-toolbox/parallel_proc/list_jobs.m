function [jm,jobs] = list_jobs(sched)
% [jm,jobs] = list_jobs(sched)
%
% sched = scheduler struct
%   .type = 'jobmanager' or 'local'
%   .name = 'jm' or 'jm2'
%   .url = 'heimdall.cluster.cresis.ku.edu';
%
% This program will optionally use the global variable
%   gRadar.sched
% which will be used to fill in any missing fields in sched
% or if sched is undefined or empty.
%
% Examples:
%   list_jobs(struct('type','jobmanager','name', ...
%     'jm2','url','heimdall.cluster.cresis.ku.edu'));
%
%   list_jobs();
%     --> Uses gRadar.sched structure
%
%   list_jobs(struct('name','jm2'));
%     --> List all jobs using gRadar.sched structure except
%         override with sched.name = 'jm2'
%
% See also list_jobs.m, get_job.m, findJob.m, findTask.m, destroy_jobs.m
 
% Get the global structure
if ~isempty(whos('global','gRadar'))
  % gRadar exists, so let us declare and use it
  global gRadar
  if ~isempty(gRadar) && isstruct(gRadar) && isfield(gRadar,'sched')
    lRadar.sched = gRadar.sched;
  end
end

if ~exist('sched','var') || isempty(sched) || ~isstruct(sched)
  if exist('lRadar','var')
    sched = lRadar.sched;
  else
    error('sched and gRadar.sched are not set');
  end
end

if ~isfield(sched,'ver')
  sched.ver = 1;
end

if ~isfield(sched,'type')
  if exist('lRadar','var') && isfield(lRadar.sched,'type')
    sched.type = lRadar.sched.type;
  else
    error('gRadar.sched.type and sched.type undefined');
  end
end
if strcmpi(sched.type,'jobmanager') && ~isfield(sched,'url')
  if exist('lRadar','var') && isfield(lRadar.sched,'url')
    sched.url = lRadar.sched.url;
  else
    error('gRadar.sched.url and sched.url undefined');
  end
end
if strcmpi(sched.type,'jobmanager') && ~isfield(sched,'name')
  if exist('lRadar','var') && isfield(lRadar.sched,'name')
    sched.name = lRadar.sched.name;
  else
    error('gRadar.sched.name and sched.name undefined');
  end
end

if strcmp(sched.type,'jobmanager')
  jm = findResource('scheduler','type',sched.type,'LookupUrl', ...
    sched.url,'name',sched.name);
elseif strcmp(sched.type,'torque')
  if isempty(sched.url)
    jm = findResource('scheduler','type',sched.type);
  else
    jm = findResource('scheduler','type',sched.type,'LookupUrl',sched.url);
  end
elseif strcmpi(sched.type,'local')
  if sched.ver == 1
    jm = findResource('scheduler','type',sched.type);
  else
    jm = parcluster();
  end
else
  error('Unsupported scheduler type %s\n', sched.type);
end

if isempty(jm)
  error('Could not find job manager');
end

jobs = get(jm,'Jobs');

if isempty(jobs)
  fprintf('No jobs\n');
else
  fprintf('%8s %12s %12s %12s %8s\n', 'Index', 'Username', 'State', ...
    'Name', 'ID');
  for ind = 1:length(jobs)
    job = jobs(ind);
  
    if sched.ver == 1
      fprintf('%8d %12s %12s %12s %8d\n', ind, get(job,'UserName'), ...
        get(job,'State'), get(job,'Name'), get(job,'ID'));
    else
      fprintf('%8d %12s %12s %12s %8d\n', ind, get(job,'Username'), ...
        get(job,'State'), get(job,'Name'), get(job,'ID'));
    end
  end
end

return;

