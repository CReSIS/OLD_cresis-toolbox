function run_sar(varargin)

% function sim.run_sar
%
% Script for running sar.m on FullSim (usually just used for debugging).
%
% Authors: John Paden, Hara Madhav Talasila
%
% See also: run_master.m, master.m, run_sar.m, sar.m, sar_task.m,
%   sar_coord_task.m

param_fn = [];

switch nargin
  case 1
    run_en = varargin{1};
  case 2
    run_en = varargin{1};
    param_fn = varargin{2};
  otherwise
    run_en = 0;
end

%% User Setup
% =====================================================================
if isempty(param_fn)
  % param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/snow/2012_Greenland_P3sim/20120330/param.mat';
  % param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/snow/2013_Greenland_P3sim/20130327/param.mat';
  % param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/snow/2018_Antarctica_DC8sim/20181010/param.mat';
  % param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/rds/2014_Greenland_P3sim/20140325/param.mat';
  % param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/rds/2018_Greenland_P3sim/20180429/param.mat';
  param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/rds/2014_Greenland_P3sim/20140410/param.mat';
  param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/rds/2014_Greenland_P3sim/20140502/param.mat';
end

% Load parameters from the mat file
load(param_fn);

%  Overrides
param_override = [];
param_override.sar.imgs           = param.sim.imgs;
if 0
  param_override.sar.sigma_x        = 5; %2*param.sar.sigma_x;
  param_override.sar.surf_filt_dist = 50; % default 3000m
  param_override.sar.chunk_len      = 50; % default 2500m
else
  param_override.sar.sigma_x        = param.sar.sigma_x;
  param_override.sar.surf_filt_dist = 3000; % default 3000m
  param_override.sar.chunk_len      = 2500; % default 2500m
end

param_override.sar.wf_adc_pair_task_group_method = 'board'; % default 'img'
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
  
  [sar_data, hdr] = sar_load(param);
  
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
  h_axes_slice = [];
  h_cb = [];
  N_imgs = length(param.sar_load.imgs);
  
  % choose slice index other than global max
  if 1
    slice_idx_override = 0;
  else
    slice_idx_override = 1;
    slice_idx_row_override = [n n n];
    slice_idx_col_override = [n n n];
  end
  leg_str_slice_row = [];
  leg_str_slice_col = [];
  
  for img = 1:N_imgs % support only one for now
    for wf_adc = 1:size(param.load_data.imgs{img},1)
      wf = param.load_data.imgs{img}(wf_adc,1);
      adc = param.load_data.imgs{img}(wf_adc,2);
      fprintf('(%s) (wf-adc) (%d-%d)\n', datestr(now), wf, adc);
      
      YLim_min = min([YLim_min; hdr.wfs(wf).time], [], 'all');
      YLim_max = max([YLim_max; hdr.wfs(wf).time], [], 'all');
      
      if size(sar_data{img},3) > 1
        single_adc_mode = 1;
        data{img} = sar_data{img}(:,:,wf_adc); % load single adc
      else
        single_adc_mode = 0;
        data{img} = sar_data{img};
      end
      
      tmp = 20*log10(abs(data{img}));
      if lat_plot_en
        x = hdr.fcs{img}{wf_adc}.lat;
        x = hdr.lat;
      else
        x = 1:length(hdr.fcs{img}{wf_adc}.lat);
        x = 1:length(hdr.lat);
      end
      
      if single_adc_mode
        fig_h = figure('Name',fig_title);
      else
        figure(fig_h);
      end
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
      
      leg_str = [];
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
        leg_str = [leg_str {'Processed', 'Simulated'} ];
      end
      
      % Global max
      [g_max, g_max_idxs] = max(tmp(:));
      [g_max_row, g_max_col] = ind2sub(size(tmp), g_max_idxs);
      plot(x(g_max_col), hdr.wfs(wf).time(g_max_row)/1e-6,'rs');
      leg_str = [leg_str {sprintf('Max %.2f',g_max)} ];
      
      legend(leg_str);
      
      if single_adc_mode
        % linkaxes(h_axes);
        % untested: linkaxes slows down due to
        % parent bug of having a mixed combo of multiple wf-adc pairs
        % solution: untangle the plots based on the images (wf, adc)
        zoom on;
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
      
      % Slices
      if 0 || sum(cellfun(@numel,param.sim.imgs))/2 > 4
        slice_en = 0;
        continue;
      else
        slice_en = 1;
      end
      
      if slice_en
        if slice_idx_override
          slice_idx_row = slice_idx_row_override(img);
          slice_idx_col = slice_idx_col_override(img);
        else
          slice_idx_row = g_max_row;
          slice_idx_col = g_max_col;
        end
        leg_str_slice_row = [leg_str_slice_row {sprintf('%d (%d-%d) %d', img, wf, adc, slice_idx_row)} ];
        leg_str_slice_col = [leg_str_slice_col {sprintf('%d (%d-%d) %d', img, wf, adc, slice_idx_col)} ];
        
        figure(fig_h_slice);
        % slow_time (along track) slice
        h_axes_slice(1) = subplot(221);
        plot(x, tmp(slice_idx_row,:) );
        grid on; hold on; axis tight;
        if lat_plot_en
          xlabel('Latitude');
        else
          xlabel('Along-track position, rlines');
        end
        ylabel('Magnitude, dB');
        legend(leg_str_slice_row);
        title('Slice slow-time / along-track [ img (wf-adc) col idx ]');
        
        
        h_axes_slice(3) = subplot(223);
        plot(x, angle(data{img}(slice_idx_row,:)) );
        grid on; hold on; axis tight;
        if lat_plot_en
          xlabel('Latitude');
        else
          xlabel('Along-track position, rlines');
        end
        ylabel('Phase, radians');
        legend(leg_str_slice_row);
        
        % fast_time (rline) slice
        h_axes_slice(2) = subplot(222);
        plot(hdr.wfs(wf).time/1e-6, tmp(:,slice_idx_col) );
        grid on; hold on; axis tight;
        xlabel('Fast-time, us');
        ylabel('Magnitude, dB');
        title(sprintf('[wf %02d adc %02d]',wf,adc));
        legend(leg_str_slice_col);
        title('Slice fast-time / rline [ img (wf-adc) col idx ]');
        
        h_axes_slice(4) = subplot(224);
        plot(hdr.wfs(wf).time/1e-6, angle(data{img}(:,slice_idx_col)) );
        grid on; hold on; axis tight;
        xlabel('Fast-time, us');
        ylabel('Phase, radians');
        legend(leg_str_slice_col);
      end %slice_en
      
    end % wf_adc
  end % img
  
  if ~single_adc_mode
    figure(fig_h);
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
  
  if slice_en
    figure(fig_h_slice);
    linkaxes(h_axes_slice([1,3]),'x'); zoom on;
    linkaxes(h_axes_slice([2,4]),'x'); zoom on;
    try
      sgtitle(fig_title,'FontWeight','Bold','FontSize',14,'Interpreter','None');
    end
    set(findobj(fig_h_slice,'type','axes'),'FontWeight', 'Bold', 'FontSize',14);
    set(fig_h_slice, 'Position', get(0, 'Screensize'));
    %   print(gcf, '-dpng', fig_title, '-r300');
  end
  
  
end