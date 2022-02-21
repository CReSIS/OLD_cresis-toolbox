% script sim.run_sar
%
% Script for running sar.m on FullSim (usually just used for debugging).
%
% Authors: John Paden, Hara Madhav Talasila
%
% See also: run_master.m, master.m, run_sar.m, sar.m, sar_task.m,
%   sar_coord_task.m

try; hara; end;

%% User Setup
% =====================================================================

% param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/snow/2012_Greenland_P3sim/20120330/param.mat';
% param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/snow/2013_Greenland_P3sim/20130327/param.mat';
% param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/snow/2018_Antarctica_DC8sim/20181010/param.mat';
% param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/rds/2014_Greenland_P3sim/20140325/param.mat';
% param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/rds/2018_Greenland_P3sim/20180429/param.mat';
param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/rds/2014_Greenland_P3sim/20140410/param.mat';

% Load parameters from the mat file
load(param_fn);

%  Overrides
param_override = [];
param_override.sar.imgs           = param.sim.imgs;
param_override.sar.surf_filt_dist = 10; % default 3000m
param_override.sar.chunk_len      = 100; % default 2500m
param.sar.wf_adc_pair_task_group_method = 'board'; % default 'img'
% default 'img' causes error  in DDC in data_pulse_compress because if the
% simulation uses more than 1 image, then the second image in simulation
% is loaded as img=1 while wf=2 which gives incorrect values for
% wfs(wf).Nt_raw for data{img}(1:wfs(wf).Nt_raw,rlines,wf_adc)

% dbstop if error;
% param_override.cluster.type = 'torque';
% param_override.cluster.type = 'matlab';
param_override.cluster.type = 'debug';
% param_override.cluster.type = 'slurm';
% param_override.cluster.rerun_only = true;
% param_override.cluster.desired_time_per_job  = 240*60;
% param_override.cluster.cpu_time_mult  = 2;
% param_override.cluster.mem_mult  = 2;

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

%% ##################################### RUN and then LOAD
%% Run SAR
% =====================================================================

% Process the segment
ctrl_chain = {};
ctrl_chain{end+1} = sar(param,param_override);
cluster_print_chain(ctrl_chain);
[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);

% Commands to load and run from any computer (chain file contains list of file dependencies):
[ctrl_chain,chain_fn] = cluster_load_chain(chain_id);
ctrl_chain = cluster_run(ctrl_chain);

%% Load SAR data
% =====================================================================

param.sar_load.imgs = param_override.sar.imgs;

[data, hdr] = sar_load(param);

call_sign = sprintf('FullSim SAR Data %s', param.day_seg);
fig_title = sprintf('%s_%s',mfilename, call_sign);
fig_h = figure('Name',fig_title);

relative_pow_en = 0;
YLim_min = Inf;
YLim_max = -Inf;
P_min = Inf;
P_max = -Inf;
h_axes = [];
h_cb = [];
N_imgs = length(param.sar_load.imgs);

for img = 1:N_imgs % support only one for now
  for wf_adc = 1:size(param.load_data.imgs{img},1)
    wf = param.load_data.imgs{img}(wf_adc,1);
    adc = param.load_data.imgs{img}(wf_adc,2);
    
    YLim_min = min([YLim_min; hdr.wfs(wf).time], [], 'all');
    YLim_max = max([YLim_max; hdr.wfs(wf).time], [], 'all');
    x = 1:length(hdr.lat);
    tmp = 20*log10(abs(data{img}));
    
    h_axes(img) = subplot(1, N_imgs, img);
    if relative_pow_en
      tmp = tmp-max(tmp(:));
      imagesc(x, hdr.wfs(wf).time/1e-6, tmp, [-30,0] );
      cb = colorbar; cb.Label.String = 'Relative Power, dB';
    else
      imagesc(x, hdr.wfs(wf).time/1e-6, tmp);
      cb = colorbar; cb.Label.String = 'Power, dB';
      P_min = min( [P_min; min(tmp(:))] );
      P_max = max( [P_max; max(tmp(:))] );
    end
    grid on; hold on; axis tight;
    % plot(x, indi{img}.Surface/1e-6, 'x', 'Color', 'g');
    xlabel('Along-track position, rlines'); ylabel('Fast-time, us');
    title(sprintf('[wf %02d adc %02d]',wf,adc));
  end
end

linkaxes(h_axes); zoom on;
set(h_axes, 'YLim', [YLim_min YLim_max]/1e-6);
if ~relative_pow_en
  set(h_axes, 'CLim', [P_min P_max]);
end
try
  sgtitle(fig_title,'FontWeight','Bold','FontSize',14,'Interpreter','None');
end
set(findobj(fig_h,'type','axes'),'FontWeight', 'Bold', 'FontSize',14);
set(fig_h, 'Position', get(0, 'Screensize'));
%   print(gcf, '-dpng', fig_title, '-r300');
