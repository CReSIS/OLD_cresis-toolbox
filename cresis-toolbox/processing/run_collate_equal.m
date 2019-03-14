% script run_collate_equal
%
% Script for running collate_equal which takes outputs from analysis
% waveform command and finds the time delay, phase, and amplitude
% differences between different wf_adc pairs.
%
% See "Receiver equalization" wiki page for details.
%
% Author: John Paden
%
% See also: collate_equal
param_override = [];

params = read_param_xls(ct_filename_param('rds_param_2018_Greenland_P3.xls'),'',{'analysis_equal' 'analysis'});

params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20180322_03');
params = ct_set_params(params,'cmd.generic',1,'day_seg','20180322_04');
% param_override.collate_equal.debug_plots = {'final'};
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20180420_01');
param_override.collate_equal.out_dir = 'analysis';
param_override.collate_equal.img_lists = {1,2,3};
param_override.collate_equal.wf_adc_idxs = {{[1:4,6:16]},{[1:4,6:16]},{[1:4,6:16]}};
params = ct_set_params(params,'radar.wfs(1).Tsys',[0.46 -4.66 0.14 -1.77 0 -2.63 -3.38 -69.66 -75.57 -75.45 -80.42 -80.49 -75.71 -77.69 -70.53]/1e9);
params = ct_set_params(params,'radar.wfs(2).Tsys',[0.46 -4.66 0.14 -1.77 0 -2.63 -3.38 -69.66 -75.57 -75.45 -80.42 -80.49 -75.71 -77.69 -70.53]/1e9);
params = ct_set_params(params,'radar.wfs(3).Tsys',[0.46 -4.66 0.14 -1.77 0 -2.63 -3.38 -69.66 -75.57 -75.45 -80.42 -80.49 -75.71 -77.69 -70.53]/1e9);
params = ct_set_params(params,'radar.wfs(1).chan_equal_dB',[6.8 -0.6 3 0.1 0 3.5 3.9 7 3.3 4.8 6.1 6.2 4.6 3.1 6.2]);
params = ct_set_params(params,'radar.wfs(2).chan_equal_dB',[6.8 -0.6 3 0.1 0 3.5 3.9 7 3.3 4.8 6.1 6.2 4.6 3.1 6.2]);
params = ct_set_params(params,'radar.wfs(3).chan_equal_dB',[6.8 -0.6 3 0.1 0 3.5 3.9 7 3.3 4.8 6.1 6.2 4.6 3.1 6.2]);
params = ct_set_params(params,'radar.wfs(1).chan_equal_deg',[-166.2 -142.7 177 -95.9 0 -25.9 -86.5 -27.4 128.1 41.6 -46.8 43 90.7 121.3 31.6]);
params = ct_set_params(params,'radar.wfs(2).chan_equal_deg',[-166.2 -142.7 177 -95.9 0 -25.9 -86.5 -27.4 128.1 41.6 -46.8 43 90.7 121.3 31.6]);
params = ct_set_params(params,'radar.wfs(3).chan_equal_deg',[-166.2 -142.7 177 -95.9 0 -25.9 -86.5 -27.4 128.1 41.6 -46.8 43 90.7 121.3 31.6]);

% param_override.collate_equal.img_lists = {[1 2 3]};
% param_override.collate_equal.wf_adc_idxs = {{[1:4,6:16],[1:4,6:16],[1:4,6:16]}};
% 
% param_override.collate_equal.img_lists = {[1 2],[3 4],[5 6]};
% param_override.collate_equal.wf_adc_idxs = {{[1:4,6:16],[1:4,6:16]},{[1:4,6:16],[1:4,6:16]},{[1:4,6:16],[1:4,6:16]}};

% params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20180315_09');
% param_override.collate_equal.out_dir = 'analysis';
% param_override.collate_equal.img_lists = 2;

%   params.analysis.surf.ref_wf_adc = 2;
%   params.analysis.surf.wf_adc_list = [1:4];
%     params.analysis.surf.ref_wf_adc = 2;
%   params.analysis.surf.wf_adc_list = [13:16];
%   params.radar.wfs(8).rx_paths = [8 9 10 11 1 1 2 3 4 5 6 7 12 13 14 15];
% params.radar.wfs(5).chan_equal_dB = [7.7 4.9 5.6 6.9 5.5 4.1 2.6 0.4 0 3.1 3.6 6.9 5.3 5 6.9];
% params.radar.wfs(5).chan_equal_deg = [-171.4 48.6 -57.8 -177.3 -178.2 -160.4 165.2 20.9 0 161.5 152.7 -73.4 27.9 129.2 -49.8];

param_override.collate_equal.ref = 9;
param_override.collate_equal.motion_comp_en = true;
param_override.collate_equal.cmd_idx = 1;
param_override.collate_equal.chan_eq_en = true;
param_override.collate_equal.rlines = [];
% zero_surf_bin_override = 11; % Normally not used except for internal layers

%   params.collate_equal.debug_level = 1; % <-- TYPICAL SETTINGS ARE 1, 3, and 4
%param.analysis.surf.rlines = [1000:2000]; % <-- MAY WANT TO HAND MODIFY
% param.analysis.surf.motion_comp.en = true; % <-- MAY WANT TO HAND MODIFY
% param.analysis.surf.wf_adc_list = [9:16]; % <-- MAY WANT TO HAND MODIFY
% param.analysis.surf.ref_wf_adc = 9; % <-- MAY WANT TO HAND MODIFY
% param.analysis.surf.retrack.en = false; % <-- MAY WANT TO HAND MODIFY

for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
%   collate_equal(param,param_override);
  collate_equal;
end
