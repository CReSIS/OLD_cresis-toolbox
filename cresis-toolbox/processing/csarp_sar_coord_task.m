function [success] = csarp_sar_task(param)
% [success] = csarp_sar_task(param)
%
% Creates the SAR coordinate system for the entire segment.
%
% Authors: John Paden
%
% Example:
%  See run_csarp.m for how to run this function directly.
%  Normally this function is called from master.m using the param spreadsheet.
%
% See also: run_master.m, master.m, run_csarp.m, csarp.m,
%   csarp_task.m

%% Setup processing
% =====================================================================

% Get speed of light, dielectric of ice constants
physical_constants;
wgs84 = wgs84Ellipsoid('meters');

% Get the standard radar name
[~,~,radar_name] = ct_output_dir(param.radar_name);

% Load records file
records_fn = ct_filename_support(param,'','records');
records = load(records_fn);

% SAR output directory
csarp_out_dir = ct_filename_out(param, param.csarp.out_path);

%% Collect waveform information into one structure
%  - This is used to break the frame up into chunks
% =====================================================================
if strcmpi(radar_name,'mcrds')
  wfs = load_mcrds_wfs(records.settings, param, ...
    records.param_records.records.file.adcs, param.csarp);
elseif any(strcmpi(radar_name,{'acords','hfrds','hfrds2','mcords','mcords2','mcords3','mcords4','mcords5','seaice','accum2'}))
  wfs = load_mcords_wfs(records.settings, param, ...
    records.param_records.records.file.adcs, param.csarp);
elseif any(strcmpi(radar_name,{'icards'}))% add icards---qishi
  wfs = load_icards_wfs(records.settings, param, ...
    records.param_records.records.file.adcs, param.csarp);
elseif any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5'}))
  error('Not supported');
  wfs = load_fmcw_wfs(records.settings, param, ...
    records.param_records.records.file.adcs, param.csarp);
  for wf=1:length(wfs)
    wfs(wf).time = param.csarp.time_of_full_support;
    wfs(wf).freq = 1;
  end
end

%% Create the along-track axes (used for structuring processing)
% This cannot be done in advance or asynchronously because it controls how
% the data will be broken into chunks. It usually only needs to be run one
% time though.
%
% sar_ref:
%   .orig: 3xNx matrix (3 by Nx, Nx = output positions)
%   .x: 3xNx matrix, along-track unit vector
%   .z: 3xNx matrix, up unit vector (orthgonal to x)
%     Note that y, left vector, can be computed as cross of z and x
%   .along_track: 1 by Nx vector
%   .pp: spline interpolation kernel for
% =====================================================================

Lsar = c/wfs(1).fc*(param.csarp.Lsar.agl+param.csarp.Lsar.thick/sqrt(er_ice))/(2*param.csarp.sigma_x);

trajectory_param = struct('gps_source',records.gps_source, ...
  'season_name',param.season_name,'radar_name',param.radar_name,'rx_path', 0, ...
  'tx_weights', [], 'lever_arm_fh', param.csarp.lever_arm_fh);

ref = trajectory_with_leverarm(records,trajectory_param);

along_track = geodetic_to_along_track(ref.lat,ref.lon,ref.elev,Lsar);
output_along_track = along_track(1) : param.csarp.sigma_x : along_track(end);

surf_idxs = get_equal_alongtrack_spacing_idxs(along_track,param.csarp.sigma_x);
if surf_idxs(end) ~= length(along_track)
  surf_idxs(end+1) = length(along_track);
end

surf = sgolayfilt(records.surface(surf_idxs),3,round(param.csarp.surf_filt_dist / median(diff(along_track(surf_idxs)))/2)*2+1);

SAR_coord_param.type = param.csarp.mocomp.type;
SAR_coord_param.squint = [0 0 -1].';
SAR_coord_param.Lsar = Lsar;
fcs = SAR_coord_system(SAR_coord_param,ref,ref,along_track,output_along_track);

sar = [];
sar.version = 1.0;
sar.Lsar = Lsar;
sar.gps_source = records.gps_source;
sar.type = param.csarp.mocomp.type;
sar.sigma_x = param.csarp.sigma_x;
sar.along_track = along_track;
sar.surf_pp = spline(along_track(surf_idxs),surf);
sar.origin = fcs.origin;
sar.x = fcs.x;
sar.z = fcs.z;
sar.roll = fcs.roll;
sar.pitch = fcs.pitch;
sar.heading = fcs.heading;
sar.gps_time = fcs.gps_time;

sar_fn = fullfile(csarp_out_dir,'sar_coord.mat');
sar_fn_dir = fileparts(sar_fn);
if ~exist(sar_fn_dir,'dir')
  mkdir(sar_fn_dir);
end
fprintf('Saving SAR coord %s (%s)\n', sar_fn, datestr(now));
save(sar_fn,'-v6','-struct','sar','version','Lsar','gps_source','type','sigma_x','surf_pp','along_track','origin','x','z','roll','pitch','heading','gps_time');

fprintf('%s done %s\n', mfilename, datestr(now));

success = true;
