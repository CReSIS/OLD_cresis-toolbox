% script runOpsGetCrossoversWithinPolygon.m
%
% Example script for running opsGetCrossoversWithinPolygon.m.
%
% Authors: Reece Mathews

% SET THE SYSTEM
sys = 'accum';

clear param;
param.properties.location = 'arctic';
param.properties.bound = 'POLYGON((-50.076388231867725 65.2897267129677,-50.21461380890868 65.24096089671465,-50.10598573426196 65.18775402828334,-49.98510940956086 65.18672753565116,-49.930472461209064 65.23996185833009,-49.9977719704356 65.27278358141487,-50.076388231867725 65.2897267129677))';
param.properties.seasons = {'2014_Greenland_P3', '2009_Greenland_TO'};


% AUTOMATED SECTION BELOW

start = tic;
fprintf('Running opsGetCrossoversWithinPolygon ...\n');
[status,data] = opsGetCrossoversWithinPolygon(sys,param);
stop = toc(start);
fprintf('   Completed in %f seconds\n',stop);

status
data
