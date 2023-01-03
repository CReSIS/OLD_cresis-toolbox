%
% A WRAPPER FOR OPS_GET_CROSSOVER_ERROR.m
%
% YOU CAN TAKE THIS OUTPUT AND RUN OPS_PARSE_CROSSOVERS.m
%

% SET AN OUTPUT DIRECTORY FOR MAT FILE/S
output_dir = 'C:\users\mathe\Documents\Code\Matlab\CReSIS\crossovers';
%output_dir = '/cresis/scratch2/paden/mdce_tmp/crossovers/';

% SET THE SYSTEM
sys = 'rds';

% SET A PARAM FOR EACH LOCATION SET YOU WANT TO PROCESS (START AT LAST # IN PARAMS(#))
params = []; param_idx = 0;

param_idx = param_idx + 1;
params(param_idx).properties.location = 'arctic';
params(param_idx).properties.lyr_name = {'surface','bottom'};
params(param_idx).properties.season = {'2008_Greenland_TO_wise','2009_Greenland_TO_wise','2010_Greenland_DC8','2010_Greenland_TO_wise','2011_Greenland_P3',...
  '2012_Greenland_P3','2013_Greenland_P3'};

param_idx = param_idx + 1;
params(param_idx).properties.location = 'antarctic';
params(param_idx).properties.lyr_name = {'surface','bottom'};
params(param_idx).properties.season = {'2002_Antarctica_P3chile','2004_Antarctica_P3chile','2009_Antarctica_DC8','2009_Antarctica_TO','2010_Antarctica_DC8',...
  '2011_Antarctica_DC8','2011_Antarctica_TO','2012_Antarctica_DC8'};

% AUTOMATED SECTION BELOW
for pidx = 1:length(params)
  
  outputFn = fullfile(output_dir,strcat(sys,params(pidx).properties.location,'Crossovers.mat'));
  fprintf('Current parameters for crossovers:\n');
  sys
  params(pidx).properties
  start = tic;
  fprintf('Running crossovers ...\n');
  [status,data] = opsGetCrossoverError(sys,params(pidx));
  stop = toc(start);
  fprintf('   Saving %s ...\n',outputFn);
  save(outputFn,'data');
  fprintf('   Completed in %f minutes\n',stop/60);  
  
end