function [job,tasks] = get_job(sched,ID)
% [jobs,tasks] = get_job(sched,ID)
%
% sched = scheduler struct OR a jobmanager handle
%   If a struct:
%   .type = 'jobmanager' or 'local'
%   .name = 'jm' or 'jm2'
%   .url = 'heimdall.cluster.cresis.ku.edu';
%   If a jobmanager handle, should be returned from findResource (e.g.
%   list_jobs does this).
% ID = scalar integer ID
%
% This program will optionally use the global variable
%   gRadar.sched
% which will be used to fill in any missing fields in sched
% or if sched is undefined or empty.
%
% Examples:
%   ID = 776;
%   job = get_job(struct('type','jobmanager','name', ...
%     'jm2','url','heimdall.cluster.cresis.ku.edu'), ID);
%
%   get_job([], ID);
%     --> Uses gRadar.sched structure
%
%   get_job(struct('name','jm2'), ID);
%     --> List all jobs using gRadar.sched structure except
%         override with sched.name = 'jm2'
%
% See also list_jobs.m, get_job.m, findJob.m, findTask.m, destroy_jobs.m

% Get Matlab Version Information
[~,version_date] = version;
version_date = datenum(version_date);
version_two_date = datenum('Sept 15, 2014'); % Estimate of when Matlab changed versions

if isstruct(sched) || isempty(sched)
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
  elseif strcmpi(sched.type,'local')
    if (version_date < version_two_date)
      jm = findResource('scheduler','type',sched.type);
    else
      jm = parcluster();
    end
  else
    error('Unsupported scheduler type %s\n', sched.type);
  end
else
  jm = sched;
end

if isempty(jm)
  error('Could not find job manager');
end

job = jm.findJob('ID',ID);

if isempty(job)
  fprintf('Job not found, try running list_jobs or findJob\n');
  job = [];
  tasks = [];
else
  job = job(1);
  get(job)
  if (version_date < version_two_date)
    tasks = get(job,'tasks')
  else
    tasks = get(job,'Tasks')
  end
end

return;

