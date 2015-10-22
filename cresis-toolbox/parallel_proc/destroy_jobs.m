function destroy_jobs(match_param,sched)
% destroy_jobs(match_param,sched)
%
% match_param = string containing case insensitive regular expression
%   user name.  See "help regexpi". Leaving empty just lists the jobs.
%     -- OR --
%   vector of job ids to destroy
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
%   destroy_jobs('paden',struct('type','jobmanager','name', ...
%     'jm2','url','heimdall.cluster.cresis.ku.edu'));
%
%   destroy_jobs('^paden$');
%     --> Uses gRadar.sched structure and exact match of 'paden'
%
%   destroy_jobs('paden');
%     --> Uses gRadar.sched structure and match of 'paden'
%
%   destroy_jobs('.*');
%     --> Destroy all jobs using gRadar.sched structure
%
%   destroy_jobs('paden',struct('name','jm2'));
%     --> Destroy all jobs using gRadar.sched structure except
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

job_struct = [];
for ind = 1:length(jobs)
  if sched.ver == 1
    job_struct.UserName = get(jobs(ind),'UserName');
  else
    job_struct.UserName = get(jobs(ind),'Username');
  end
  job_struct.State = get(jobs(ind),'State');
  
  if exist('match_param','var') ...
      && ((ischar(match_param) && ~isempty(regexpi(job_struct.UserName,match_param))) ...
      || any(get(jobs(ind),'ID') == match_param))
    fprintf('Destroying job (%s,%s,%d)\n', job_struct.UserName, ...
      job_struct.State, get(jobs(ind),'ID'));
    destroy(jobs(ind));
  else
    fprintf('Skipping job (%s,%s,%d)\n', job_struct.UserName, ...
      job_struct.State, get(jobs(ind),'ID'));
  end
end

return;

