% script run_analysis
%
% Script for running analysis
% https://ops.cresis.ku.edu/wiki/index.php/Analysis
%
% Authors: John Paden
%
% See also: master.m, run_analysis.m, analysis.m, analysis_task.m

%% User Setup
% =====================================================================
param_override = [];

% param_override.radar.wfs(1).prepulse_H.type = 'inverse_filter';
% param_override.radar.wfs(1).prepulse_H.fn = 'bottom';

% params = read_param_xls(ct_filename_param('snow_param_2018_Antarctica_DC8.xls'),'',{'analysis_noise' 'analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'),'',{'analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'),'',{'analysis_noise','analysis'});
params = read_param_xls(ct_filename_param('snow_param_2013_Greenland_P3.xls'),'',{'analysis_noise','analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2014_Greenland_P3.xls'),'',{'analysis_noise','analysis'});


% params = ct_set_params(params,'radar.wfs.DDC_valid',1);
% params = ct_set_params(params,'radar.wfs.nz_valid',[0 1 2 3]);

% params = ct_set_params(params,'analysis.cmd{2}.start_time','s = es.Tstart+10*es.dt');
% params = ct_set_params(params,'analysis.cmd{2}.stop_time','s = es.Tend-10*es.dt');

% params = ct_set_params(params,'analysis.cmd{4}.start_time','s = es.Tstart+10*es.dt');
% params = ct_set_params(params,'analysis.cmd{4}.stop_time','s = es.Tend-10*es.dt');

% param_override.analysis.cmd{1}.start_time = 'es.Tstart+10*es.dt';
% param_override.analysis.cmd{1}.stop_time = 'es.Tend-10*es.dt';
% 
% param_override.analysis.cmd{2}.start_time = 'es.Tstart+10*es.dt';
% param_override.analysis.cmd{2}.stop_time = 'es.Tend-10*es.dt';

% Example to run a specific segment and frame by overriding parameter spreadsheet values
params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20120314_05');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20120330_04');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20130420_01');
params = ct_set_params(params,'cmd.generic',1,'day_seg','20130328_01'); %random test for reuse_debug in collate_coh_noise
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20140428_01');

% params = ct_set_params(params,'analysis.cmd{2}.num_sam_hint',{11}); % for cluster mem and time calculations

% params = ct_set_params(params,'analysis.cmd{2}.pulse_comp',1); % done in the sheet
% cmd.pulse_comp

dbstop if error;
param_override.cluster.type = 'torque';
% param_override.cluster.type = 'matlab';
% param_override.cluster.type = 'debug';
% param_override.cluster.rerun_only = true;
% param_override.cluster.desired_time_per_job  = 240*60;
param_override.cluster.cpu_time_mult  = 10;
param_override.cluster.mem_mult  = 5;

% %since we know current ftg test has 3 waveforms
% param_override.radar.wfs(1).ftg_comp_en = 0; % Set non-empty or no field for ftg compensation
% param_override.radar.wfs(2).ftg_comp_en = 0; % Set non-empty or no field for ftg compensation
% param_override.radar.wfs(3).ftg_comp_en = 0; % Set non-empty or no field for ftg compensation

% param_override.radar.wfs(1).TTL_start = 424;
% param_override.radar.wfs(1).TTL_length = 239;
% param_override.radar.wfs(2).TTL_start = 424;
% param_override.radar.wfs(2).TTL_length = 350;
% param_override.radar.wfs(3).TTL_start = 424;
% param_override.radar.wfs(3).TTL_length = 739;
% 
% param_override.radar.wfs(1).record_start = 1328;
% param_override.radar.wfs(1).record_stop = 5776;
% param_override.radar.wfs(2).record_start = 1328;
% param_override.radar.wfs(2).record_stop = 5776;
% param_override.radar.wfs(3).record_start = 1328;
% param_override.radar.wfs(3).record_stop = 5776;
% 
% param_override.radar.TTL_clock = 1e9/18;
% 
% param_override.radar.ftg_dir = 'X:\ct_data\ct_tmp\waveform\rds\2019_Greenland_P3\gain_correction_dir';

%% Automated Section
% =====================================================================

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

% Process each of the segments
ctrl_chain = {};
for param_idx = 1:length(params)
  param = params(param_idx);
  if isfield(param.cmd,'generic') && ~iscell(param.cmd.generic) && ~ischar(param.cmd.generic) && param.cmd.generic
    ctrl_chain{end+1} = analysis(param,param_override);
  end
end

cluster_print_chain(ctrl_chain);

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);

