function run_sar(varargin)

% function sim.run_sar
%
% Script for running sar.m on FullSim (usually just used for debugging).
%
% Authors: John Paden, Hara Madhav Talasila
%
% See also: run_master.m, master.m, run_sar.m, sar.m, sar_task.m,
%   sar_coord_task.m

switch nargin
  case 1
  run_en = varargin{1};
  otherwise
  run_en = 0;
end

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

param_override.sar.start_eps = 1; % for now
% for later: follow John's email


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

if run_en %% ##################################### RUN and then LOAD
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
  
else
  %% Load SAR data
  % =====================================================================
  
  [c, WGS84] = physical_constants('c', 'WGS84');
  
  param.sar_load.imgs = param_override.sar.imgs;
  
  [data, hdr] = sar_load(param);
  
  call_sign = sprintf('FullSim SAR Data %s', param.day_seg);
  fig_title = sprintf('%s_%s',mfilename, call_sign);
  fig_h = figure('Name',fig_title);
  fig_h_slice = figure('Name',fig_title);
  
  relative_pow_en = 0;  % #############################
  
  if any( strcmpi(param.target.type, {'point', 'points'}) )
    point_target_marker_en = 1;
  else
    point_target_marker_en = 0;
  end
  
  if point_target_marker_en && param.sim.north_along_track_en
    lat_plot_en = 1;
  else
    lat_plot_en = 0;
  end
  
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
      tmp = 20*log10(abs(data{img}));
      if lat_plot_en
        x = hdr.fcs{img}{wf_adc}.lat;
        x = hdr.lat;
      else
        x = 1:length(hdr.fcs{img}{wf_adc}.lat);
        x = 1:length(hdr.lat);
      end
      
      figure(fig_h);
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
      if lat_plot_en
        xlabel('Latitude');
      else
        xlabel('Along-track position, rlines');
      end
      ylabel('Fast-time, us');
      title(sprintf('[wf %02d adc %02d]',wf,adc));
      
      if point_target_marker_en
        e_lon = abs(bsxfun(@minus, hdr.fcs{img}{wf_adc}.lon, param.target.lon'));
        e_lat = abs(bsxfun(@minus, hdr.fcs{img}{wf_adc}.lat, param.target.lat'));
        [dev_lon, idx_lon] = min(e_lon,[],2);
        [dev_lat, idx_lat] = min(e_lat,[],2);
        % e_lon, dev_lon, idx_lon are not relevant for northward flightpath
        % use idx_lat for marker plots
        range_est = distance_geodetic(hdr.fcs{img}{wf_adc}.lat(idx_lat), ...
          hdr.fcs{img}{wf_adc}.lon(idx_lat), hdr.fcs{img}{wf_adc}.elev(idx_lat), ...
          param.target.lat, param.target.lon, param.target.elev, WGS84.ellipsoid) ;
        TWTT_est = range_est * 2/c;
        plot(hdr.fcs{img}{wf_adc}.lat(idx_lat),TWTT_est/1e-6,'kx');
        plot(param.target.lat, TWTT_est/1e-6,'ko');
        legend({'Loaded', 'Simulated'});
      end
      
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
  
end