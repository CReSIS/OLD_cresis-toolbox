% =========================================================================
% OPS COMMAND SCRIPT
%
% This script is used throughout the OPS MATLAB API and sets up some basic
% parameters and variables used throughout the script.
%
% Please set the user input below.
%
% Authors: Kyle W. Purdon, Trey Stafford
%
% =========================================================================

%% USER INPUT

global gRadar;

% THE OPS PROFILER WILL RUN AND RETURN PROFILING LOGS
if isfield(gRadar,'ops') && isfield(gRadar.ops,'profileCmd') && ~isempty(gRadar.ops.profileCmd)
  gOps.profileCmd = gRadar.ops.profile_cmd;
else
  gOps.profileCmd = false;
end

if isfield(gRadar,'ops') && isfield(gRadar.ops,'url') && ~isempty(gRadar.ops.url)
  gOps.sysUrl = gRadar.ops.url;
else
  gOps.sysUrl = 'https://ops.cresis.ku.edu/'; % Read-only for outside of CReSIS
  % gOps.sysUrl = 'http://ops.cresis.ku.edu/'; % Use from within CReSIS
  % gOps.sysUrl = 'http://192.168.111.222/'; % Virtual box setup
end

%% AUTOMATED SECTION (DO NOT MODIFY)

gOps.dbUser = '';
gOps.dbPswd = '';

gOps.serverUrl = strcat(gOps.sysUrl,'ops/');
gOps.geoServerUrl = strcat(gOps.sysUrl,'geoserver/');

if gOps.profileCmd
    web(strcat(gOps.serverUrl(1:end-4),'profile-logs/'));
end
