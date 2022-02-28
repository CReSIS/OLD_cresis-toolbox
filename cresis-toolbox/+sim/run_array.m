% script sim.run_array
%
% Script for running array.m on FullSim(usually just used for debugging).
%
% Authors: John Paden, Hara Madhav Talasila
%
% See also: run_master.m, master.m, run_array.m, array.m, load_sar_data.m,
% array_proc.m, array_task.m, array_combine_task.m

try; hm; end;

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
param_override.array.imgs           = param.sim.imgs;
% param_override.sar.surf_filt_dist = 10; % default 3000m
param_override.array.chunk_len      = 100; % default 2500m

param_override.sar.imgs             = param.sim.imgs;

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

if 0 %% ##################################### RUN and then LOAD
  %% Run ARRAY
  % =====================================================================
  
  % Process the segment
  ctrl_chain = {};
  ctrl_chain{end+1} = array(param,param_override);
  cluster_print_chain(ctrl_chain);
  [chain_fn,chain_id] = cluster_save_chain(ctrl_chain);
  
  % Commands to load and run from any computer (chain file contains list of file dependencies):
  [ctrl_chain,chain_fn] = cluster_load_chain(chain_id);
  ctrl_chain = cluster_run(ctrl_chain);
  
else
  %% Load ARRAY data
  % =====================================================================
  
  [c, WGS84] = physical_constants('c', 'WGS84');  
  
  if isfield(param.sim,'frame_idx')
    frm = param.sim.frm;
  else
    frm = 1;
  end
  
  out_dir = ct_filename_out(param, 'standard');
  
  relative_pow_en = 0; % #############################
  
  if strcmpi(param.target.type, 'point')
    point_target_marker_en = 1;
  else
    point_target_marker_en = 0;
  end
  
  if point_target_marker_en && param.sim.north_along_track_en
    lat_plot_en = 1;
  else
    lat_plot_en = 0;
  end
  
  YLim_min = 0;
  YLim_max = 0;
  h_axes = 0;
  N_imgs = length(param_override.array.imgs); %length(param.sim.imgs);
  
  for idx = 1:N_imgs+1
    img = idx-1;
    
    if img == 0 % Combined image
      out_fn = fullfile(out_dir, sprintf('Data_%s_%03d.mat', param.day_seg, frm));
      full = load(out_fn);
      
      YLim_min = min([YLim_min; full.Time], [], 'all');
      YLim_max = max([YLim_max; full.Time], [], 'all');
      tmp = 20*log10(abs(full.Data));
      if lat_plot_en
        x = full.Latitude;
      else
        x = 1:length(full.Latitude);
      end
      
      call_sign = sprintf('FullSim Array Data %s', param.day_seg);
      fig_title = sprintf('%s_%s',mfilename, call_sign);
      fig_h = figure('Name',fig_title);
      
      h_axes(idx) = subplot(1, N_imgs+1, idx);
      if relative_pow_en
        tmp = tmp-max(tmp(:));
        imagesc(x, full.Time/1e-6, tmp, [-30,0] );
        cb = colorbar; cb.Label.String = 'Relative Power, dB';
      else
        imagesc(x, full.Time/1e-6, tmp);
        cb = colorbar; cb.Label.String = 'Power, dB';
      end
      grid on; hold on; axis tight;
      if lat_plot_en
        xlabel('Latitude');
      else
        xlabel('Along-track position, rlines');
      end
      ylabel('Fast-time, us');
      title('Combined');
      
      if point_target_marker_en
        e_lon = abs(full.Longitude - param.target.lon);
        e_lat = abs(full.Latitude - param.target.lat);
        [dev_lon, idx_lon] = min(e_lon);
        [dev_lat, idx_lat] = min(e_lat);
        % e_lon, dev_lon, idx_lon are not relevant for northward flightpath
        % use idx_lat for marker plots
        range_est = distance_geodetic(full.Latitude(idx_lat), ...
          full.Longitude(idx_lat), full.Elevation(idx_lat), ...
          param.target.lat, param.target.lon, param.target.elev, WGS84.ellipsoid) ;
        TWTT_est = range_est * 2/c;
        plot(full.Latitude(idx_lat),TWTT_est/1e-6,'kx');
        plot(param.target.lat, TWTT_est/1e-6,'ko');
      end
      
    else % individual images
      wf_adc = 1;
      wf = param.load_data.imgs{img}(wf_adc,1);
      adc = param.load_data.imgs{img}(wf_adc,2);
      
      out_fn = fullfile(out_dir, sprintf('Data_img_%02d_%s_%03d.mat', img, param.day_seg, frm));
      indi{img} = load(out_fn);
      
      YLim_min = min([YLim_min; indi{img}.Time], [], 'all');
      YLim_max = max([YLim_max; indi{img}.Time], [], 'all');
      if lat_plot_en
        x = indi{img}.Latitude;
      else
        x = 1:length(indi{img}.Latitude);
      end
      tmp = 20*log10(abs(indi{img}.Data));
      
      h_axes(idx) = subplot(1, N_imgs+1, idx);
      if relative_pow_en
        tmp = tmp-max(tmp(:));
        imagesc(x, indi{img}.Time/1e-6, tmp, [-30,0] );
        cb = colorbar; cb.Label.String = 'Relative Power, dB';
      else
        imagesc(x, indi{img}.Time/1e-6, tmp);
        cb = colorbar; cb.Label.String = 'Power, dB';
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
        e_lon = abs(indi{img}.Longitude - param.target.lon);
        e_lat = abs(indi{img}.Latitude - param.target.lat);
        [dev_lon, idx_lon] = min(e_lon);
        [dev_lat, idx_lat] = min(e_lat);
        % e_lon, dev_lon, idx_lon are not relevant for northward flightpath
        % use idx_lat for marker plots
        range_est = distance_geodetic(indi{img}.Latitude(idx_lat), ...
          indi{img}.Longitude(idx_lat), indi{img}.Elevation(idx_lat), ...
          param.target.lat, param.target.lon, param.target.elev, WGS84.ellipsoid) ;
        TWTT_est = range_est * 2/c;
        plot(indi{img}.Latitude(idx_lat),TWTT_est/1e-6,'kx');
        plot(param.target.lat, TWTT_est/1e-6,'ko');
        legend({'Loaded', 'Simulated'});
      end
      
    end
    
  end
  
  linkaxes(h_axes); zoom on;
  set(h_axes, 'YLim', [YLim_min YLim_max]/1e-6);
  try
    sgtitle(fig_title,'FontWeight','Bold','FontSize',14,'Interpreter','None');
  end
  set(findobj(fig_h,'type','axes'),'FontWeight', 'Bold', 'FontSize',14);
  set(fig_h, 'Position', get(0, 'Screensize'));
  %   print(gcf, '-dpng', fig_title, '-r300');
  
end