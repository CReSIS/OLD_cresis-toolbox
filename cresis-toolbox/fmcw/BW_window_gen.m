function BW_window_gen(param,param_override)
% BW_window_gen(param,param_override)
%
% BW_window_gen DESCRIPTION
%
% param: struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override: parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Example:
%  See run_BW_window_gen.m for how to run this function directly.
%  This function may be called from the run_master.m script using the
%  param spreadsheet and the cmd.generic column.
%
% Authors: FILL_IN_AS_NEEDED
%
% See also: run_BW_window_gen.m, BW_window_gen.m, BW_window_gen_task.m, FILL_IN_AS_NEEDED

%% General Setup
% =====================================================================

param = merge_structs(param, param_override);

% fprintf('=====================================================================\n');
% fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
% fprintf('=====================================================================\n');

%% Input checks
% =====================================================================

% Load records
records = records_load(param);

% Check BW_window_gen parameters
if ~isfield(param,'BW_window_gen') || isempty(param.BW_window_gen)
  param.BW_window_gen = [];
end

if ~isfield(param.BW_window_gen,'BW_guard') || isempty(param.BW_window_gen.BW_guard)
  param.BW_window_gen.BW_guard = [5 5];
  param.BW_window_gen.BW_guard_units = '%';
end

if ~isfield(param.BW_window_gen,'BW_guard_units') || isempty(param.BW_window_gen.BW_guard_units)
  param.BW_window_gen.BW_guard_units = 'Hz';
end

if ~isfield(param.BW_window_gen,'imgs') || isempty(param.BW_window_gen.imgs)
  param.BW_window_gen.imgs = {[1 1]};
end
param.load.imgs = param.BW_window_gen.imgs;

% stretch_deramp: LO is extended so that RF is fully supported
if ~isfield(param.BW_window_gen,'stretch_deramp') || isempty(param.BW_window_gen.stretch_deramp)
  param.BW_window_gen.stretch_deramp = false;
end

%% Setup processing
% =====================================================================
param.radar.wfs = data_load_wfs(param,records);

wf = 1;
BW = abs(param.radar.wfs(wf).f1-param.radar.wfs(wf).f0);
BW_window_ideal = sort([param.radar.wfs(wf).f0 param.radar.wfs(wf).f1]);
fs_dec = param.radar.fs / lcm_vector(param.radar.wfs(wf).DDC_valid);
dt = 1/fs_dec;
% 2* is to ensure even number of samples
df = 2 * abs(param.radar.wfs(wf).chirp_rate) * dt;

% Add BW guards/margins (due to roll off on the edges of the spectrum
% Also remove the non-overlapped portion of the spectrum (this is where the
% LO and RF signals do not overlap)
if param.BW_window_gen.BW_guard_units(1) == '%'
  param.BW_window_gen.BW_guard = BW*param.BW_window_gen.BW_guard/100;
end
BW_window(1) = BW_window_ideal(1) + param.BW_window_gen.BW_guard(1);
BW_window(2) = BW_window_ideal(2) - ~param.BW_window_gen.stretch_deramp*(max(param.radar.wfs(wf).nz_valid)+1) * param.radar.fs/2 - param.BW_window_gen.BW_guard(2);

% Ensure that the BW window is a multiple of df
BW_window = round(BW_window(1)/df)*df + [0 round(diff(BW_window)/df)*df];

fprintf('%s\t[%.12g %.12g]\n',param.day_seg, BW_window);
