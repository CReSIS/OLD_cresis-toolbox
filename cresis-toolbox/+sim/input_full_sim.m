% Script sim.input_full_sim
%
% Generates simulated data from an extracted flightline
% Save simulated raw_data, params, records, frames
%
% Authors: John Paden, Hara Madhav Talasila
%
% See also sim.run_full_sim sim.flightline_extract, run_load_data (example 7)

try; hara; end;

param=[];
data = [];
frames = [];
records = [];

%% INPUTS to sim.flightline_extract.m
% =========================================================================
% REQUIRED
% Images to sim/load  {[ones([7 1]),(1:7).'], [2*ones([7 1]),(1:7).'], [3*ones([7 1]),(1:7).']};

% param.sim.radar_name  = 'snow';
% param.sim.season_name = '2012_Greenland_P3'; %'2018_Antarctica_DC8';
% param.sim.day_seg     = '20120330_04'; % '20181010_01';
%
% param.sim.radar_name  = 'snow';
% param.sim.season_name = '2013_Greenland_P3';
% param.sim.day_seg     = '20130327_02';
%
% param.radar_name      = 'snow';
% param.season_name     = '2013_Greenland_P3';
% param.day_seg         = '20130327_02';
%
% param.sim.radar_name  = 'snow';
% param.sim.season_name = '2016_Greenland_Polar6';
% param.sim.day_seg     = '20160414_01';
%
% param.sim.radar_name  = 'snow';
% param.sim.season_name = '2018_Antarctica_DC8';
% param.sim.day_seg     = '20181010_01';

% param.sim.radar_name  = 'rds';
% param.sim.season_name = '2014_Greenland_P3';
% param.sim.day_seg     = '20140325_05';
% param.sim.imgs        = {[1 5]};
% param.sim.imgs        = {[1 5], [2 5] , [3 5]};

param.sim.radar_name  = 'rds';
param.sim.season_name = '2014_Greenland_P3';
param.sim.day_seg     = '20140410_01';
param.sim.frame_idx   = [];
param.sim.imgs        = {[1 5]};
param.sim.imgs        = {[1 5], [2 5] , [3 5]};
% %
% param.sim.imgs        = {[ones(3,1), (4:6).'], [2*ones(3,1), (4:6).'] , [3*ones(3,1), (4:6).']};
% param.sim.imgs        = {[ones(3,1), (4:6).']};
% param.sim.imgs        = {[ones([7 1]),(2:8).'], [2*ones([7 1]),(2:8).'], [3*ones([7 1]),(2:8).']};

% param.sim.radar_name  = 'rds';
% param.sim.season_name = '2018_Greenland_P3';
% param.sim.day_seg     = '20180429_01';
% param.sim.imgs        = {[1 9]};
% param.sim.imgs        = {[1 9], [3 9], [5 9]};

% =========================================================================
% Optional

% Target offset
% x: for any distance <Lsar from the begining of flightpath
% y: positive(left to flighline) negative(left to flightline)
% z: positive(above flightline) negative(below flightline)
%
% param.target.offsets.x = []; % usually Lsar/2 or [0,Lsar]
% param.target.offsets.y = []; % usually [0] or +ve/-ve value
% param.target.offsets.z = []; % usually -ve


if 1
  % Northward flight instead of actual trajectory
  % single target in the center
  param.sim.north_along_track_en = 1;
  
  param.target.offsets.x = @(Lsar) 500;  % for any distance <Lsar
  param.target.offsets.x = @(Lsar) Lsar/2; % middle of flighpath
  param.target.offsets.y = 0;
  param.target.offsets.z = -500;
  param.target.offsets.z = -base2dec('KU', 36); % about  5 us TWTT for 750 meter
  % param.target.offsets.z = -base2dec('157', 36); % about 9.89 us TWTT for 1483 meter
  
  
elseif 0
  % Northward flight instead of actual trajectory
  % multiple targets in the center rline
  param.sim.north_along_track_en = 1;
  param.target.offsets.x = @(Lsar) [Lsar/2; Lsar/2];  % for any distance <Lsar
  param.target.offsets.y = [0; 0];
  param.target.offsets.z = -base2dec({'KU';'157';'17A'}, 36);
  % about 5 us TWTT for 750 meter, 9.89 us TWTT for 1483 meter
  
elseif 1
  % Northward flight instead of actual trajectory
  % multiple targets anywhere
  param.sim.north_along_track_en = 1;
  param.target.offsets.x = @(Lsar) Lsar .* ones(3,3,3) .* [1/4 1/2 3/4] ;  % 3 rlines along the way % repmat( ([1/4; 1/2; 3/4]*[1 1 1])' , [1,1,3] )
  param.target.offsets.y = permute( ones(3,3,3) .* [-100 0 50], [3,1,2] ); % 50m right, center and left
  param.target.offsets.z = ones(3,3,3) .* -base2dec({'HM';'TA';'17A'}, 36);
  % TWTT us for meter: about 5 us 634 m, 6.17 us 1054 m, 10.4 us 1558 m
  
end

param.sim.debug_plots_en = 1;

% =========================================================================

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.sim.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% flightline_extract

[ param, records, records2, frames, exec_good ] = sim.flightline_extract(param);

if ~exec_good
  fprintf('flightline_extract executed incompletely\n');
  return;
else
  fprintf('(%s) FullSim Flightline \t--EXTRACTED\n',  datestr(now));
  fprintf('=====================================================================\n');
end

clear exec_good;

%% Simulator

[c, WGS84] = physical_constants('c', 'WGS84');
debug_plot_en  = 0;
raw_data = [];

Nx = length(param.gps.elev); % rec_len
Ntgts = length(param.target.x);
Nimgs = sum(cellfun(@numel,param.sim.imgs))/2;
fprintf('[ %d ] images, [ %d ] targets, [ %d ] records\n', Nimgs, Ntgts, Nx);

wfs = param.radar.wfs;

[output_dir,radar_type,radar_name] = ct_output_dir(param.sim.radar_name);

if param.sim.debug_plots_en
  figure(369);
  plot(records2.lon, records2.lat, '.', 'LineWidth', 5);
  hold on; grid on;
  xlabel('Longitude, Degrees');
  ylabel('Latitude, Degrees');
  leg_str= {'Reference'};
  legend(leg_str);
  title('Trajectories');
  line_markers = ['x', 'd', 'h']; % assuming only 3 wfs max
end

% Calculate raw_data matrix cell(s)
for img = 1:length(param.sim.imgs)
  for wf_adc = 1:size(param.sim.imgs{img},1)
    
    wf = param.sim.imgs{img}(wf_adc,1);
    adc = param.sim.imgs{img}(wf_adc,2);
    rx = param.radar.wfs(wf).rx_paths(adc);
    
    % Create actual trajectory #### from data_load %% Create traj
    trajectory_param = struct('gps_source',records.gps_source, ...
      'season_name',param.season_name,'radar_name',param.radar_name, ...
      'rx_path', rx, ...
      'tx_weights', param.radar.wfs(wf).tx_weights, ...
      'lever_arm_fh', param.radar.lever_arm_fh);
    [gps,lever_arm_val] = trajectory_with_leverarm(records,trajectory_param);
    
    [gps.x, gps.y, gps.z] = geodeticD2ecef(gps.lat, gps.lon, gps.elev, WGS84.ellipsoid);
    
    if param.sim.debug_plots_en
      figure(369);
      plot(gps.lon, gps.lat, line_markers(wf), 'LineWidth', img);
      leg_str = [leg_str {sprintf('wfs %02d adc %02d',wf,adc)}];
      legend(leg_str);
    end
    
    % for phase centers
    lever_arm_val = [+1; -1; -1] .* lever_arm_val;
    xx = lever_arm_val(1);
    yy = lever_arm_val(2);
    zz = lever_arm_val(3);
    
    % figure('Name', 'Phase Centers'); hold on;
    %permute( [max(param.la.pc(1,:,:)) - min(param.la.pc(1,:,:))]/10 ,[3,2,1])
    if img == 1 &&  wf_adc == 1
      h_fig = get(groot, 'Children');
      h_fig = figure(h_fig(find(strcmpi({h_fig.Name}, 'Phase Centers'),1)));
      %set(gca,'ColorOrderIndex',1);
      marker_color = colororder; % get 7x3 (colors)x(RGB)
      marker_color = zeros(7,3);
      marker_color_N = size(marker_color,1);
    end
    
    figure(h_fig);
    marker_color_idx = rem(rx-1, 7) +1; % fixed color for a rx path
    subplot(221);
    plot3(xx,yy,zz, '+', 'Color', marker_color(marker_color_idx,:), 'LineWidth', 1.5);
    text(xx+0.0004,yy,zz-0.03, sprintf('%d',rx), 'FontSize',8, 'Color', marker_color(marker_color_idx,:) );
    subplot(222); hold on;
    plot3(xx,yy,zz, '+', 'Color', marker_color(marker_color_idx,:), 'LineWidth', 1.5);
    text(xx+0.0004,yy,zz-0.03, sprintf('%d',rx), 'FontSize',8, 'Color', marker_color(marker_color_idx,:) );
    subplot(223); hold on;
    plot3(xx,yy,zz, '+', 'Color', marker_color(marker_color_idx,:), 'LineWidth', 1.5);
    text(xx+0.0004,yy,zz-0.03, sprintf('%d',rx), 'FontSize',8, 'Color', marker_color(marker_color_idx,:) );
    subplot(224); hold on;
    plot3(xx,yy,zz, '+', 'Color', marker_color(marker_color_idx,:), 'LineWidth', 1.5);
    text(xx+0.0004,yy,zz-0.03, sprintf('%d',rx), 'FontSize',8, 'Color', marker_color(marker_color_idx,:) );
    
    
    time_raw    = wfs(wf).time_raw;
    Tpd         = wfs(wf).Tpd;
    f0          = wfs(wf).f0;
    fc          = wfs(wf).fc;
    chirp_rate  = wfs(wf).chirp_rate;
    rx_tukey    = wfs(wf).tukey;
    fs          = param.radar.fs;
    dt          = 1/fs;
    if ( time_raw(2)-time_raw(1) ) - 1/fs > 1e-18
      keyboard;
    end
    fcs = fcs_local(gps);
    theta = nan(Ntgts,Nx);
    phi   = nan(Ntgts,Nx);
    gamma = nan(Ntgts,Nx);
    
    % vect_target = [Ntgts, Nx, 3] matrix; 3rd dim is (x,y,z)
    vect_target = reshape( [ (param.target.x.' - gps.x), (param.target.y.' - gps.y), (param.target.z.' - gps.z ) ], [Ntgts, Nx, 3]);
    range = sqrt( sum(vect_target.^2, 3) );
    if 0
      % point target for now
      % range = sqrt( (gps.x - param.target.x).^2 + (gps.y - param.target.y).^2 + (gps.z - param.target.z).^2 );
      range_test = sqrt( ...
        + ( bsxfun(@minus, gps.x, param.target.x') ).^2 ...
        + ( bsxfun(@minus, gps.y, param.target.y') ).^2 ...
        + ( bsxfun(@minus, gps.z, param.target.z') ).^2 );
    end
    
    for tgt = 1:Ntgts
      % theta = Elevation angle ( -fcs.Z -> vect_target)
      % [ usually 0 (nadir) to 90 (velocity vector fcs.X or XY plane or yaw plane) ]
      % back to 180 (zenith)
      % phi   = Azimuth angle   ( +fcs.Y -> vect_target)
      % [ usually 0 (left) to 90 (nadir -fcs.Z or XZ plane or pitch plane) to 180 (right) ]
      % back to 180 left
      % gamma = Angle (+fcs.X -> vect_target)
      % [usually 0 (velocity vector or fcs.X) to 90 (nadir or YZ plane or roll plane) to 180 (aft or -fcs.X)
      % atan2d range =[-180, 180], but num>=0 because of vecnorm
      % atan2d range in this case = [0, 180];
      theta_num     = vecnorm( cross(-fcs.Z.', squeeze(vect_target(tgt,:,:)), 2), 2, 2);
      phi_num       = vecnorm( cross(+fcs.Y.', squeeze(vect_target(tgt,:,:)), 2), 2, 2);
      gamma_num     = vecnorm( cross(+fcs.X.', squeeze(vect_target(tgt,:,:)), 2), 2, 2);
      theta_den     = dot(-fcs.Z.', squeeze(vect_target(tgt,:,:)), 2);
      phi_den       = dot(+fcs.Y.', squeeze(vect_target(tgt,:,:)), 2);
      gamma_den     = dot(+fcs.X.', squeeze(vect_target(tgt,:,:)), 2);
      theta(tgt,:)  = atan2d(theta_num, theta_den);
      phi(tgt,:)    = atan2d(phi_num, phi_den);
      gamma(tgt,:)  = atan2d(gamma_num, gamma_den);
    end
    
    clear fcs vect_target theta_num phi_num theta_den phi_den
    
    % figure;imagesc(gamma); colorbar
    % figure;imagesc(theta); colorbar
    % figure;imagesc(phi); colorbar
    
    twtt = 2*range/c;
    
    Nt = wfs(wf).Nt_raw;
    raw_data{img}(:,:,wf_adc) = zeros(Nt,Nx);
    
    for tgt = 1:Ntgts
      td    = twtt(tgt,:);
      fb    = chirp_rate*td;
      time  = time_raw - td; % ################
      [cls_td, cls_time_idx] = min(td);
      cls_range = cls_td*c/2;
      
      t_zero = time_raw( find(time_raw>=0,1) );
      
      time_raw_bound_idxs = [ find(time_raw<=cls_td,1,'last') find(time_raw>=cls_td,1,'first') ];
      time_raw_bounds = time_raw(time_raw_bound_idxs);
      
      [cls_td_offset, cls_td_round_idx] = min(abs(time_raw-cls_td));
      cls_td_round = time_raw(cls_td_round_idx); %cls_td_round = dt*round(cls_td/dt);
      cls_range_round = cls_td_round*c/2;
      cls_td_err = cls_td_round - cls_td;
      cls_range_err = cls_range_round - cls_range;
      
      if 0
        figure;
        plot(time_raw,time_raw, 'x','LineWidth',2);
        hold on;
        plot(td,td, 'o','LineWidth',2);
        plot(cls_td,cls_td, 'x','LineWidth',4);
        plot(cls_td_round,cls_td_round, 'd','LineWidth',4);
      end
      
      if strcmpi(radar_type,'deramp') % deramp
        % fill in for snow
      else % pulsed
        raw_data{img}(:,:,wf_adc) = raw_data{img}(:,:,wf_adc) + tukeywin_cont(time/Tpd-0.5,rx_tukey) ...
          * 1 .* exp(1i*(2*pi*f0*time + pi*chirp_rate*time.^2 )); % Nt x Nrx;
        if 0 %%% #######################################################################################
          raw_data{img}(:,:,wf_adc) = raw_data{img}(:,:,wf_adc) .*(sqrt(4*pi*cls_range^2)) ;
        end
      end
    end % for tgt
    
    if 0 %%% #######################################################################################
      raw_data{img}(:,:,wf_adc) = raw_data{img}(:,:,wf_adc) ./(sqrt(Nt)) ./Ntgts ;
    end
    
    param.traj{img,wf_adc}      = gps;
    param.sim.range{img,wf_adc} = range;
    param.sim.twtt{img,wf_adc}  = twtt;
    
    % radiometric_corr_dB
    if strcmpi(radar_type,'deramp') % deramp
      % fill in for snow
    else % pulsed
      raw_data{img}(:,:,wf_adc) = raw_data{img}(:,:,wf_adc) * 1; %./ 10.^(wfs(wf).system_dB/20);  % Nt x Nrx
    end
    
    % Do the reverse of what happens to ~raw_data in data_load
    % Apply channel compensation, presum normalization, and constant
    % receiver gain compensation
    if strcmpi(radar_type,'deramp')
      chan_equal = 1;
    else
      chan_equal = 10.^(param.radar.wfs(wf).chan_equal_dB(param.radar.wfs(wf).rx_paths(adc))/20) ...
        .* exp(1i*param.radar.wfs(wf).chan_equal_deg(param.radar.wfs(wf).rx_paths(adc))/180*pi);
    end
    
    if length(wfs(wf).system_dB) == 1
      % Only a single number is provided for system_dB so apply it to all
      % receiver paths
      mult_factor{img}(:,:,wf_adc) = single(wfs(wf).quantization_to_V(adc) ...
        / (10.^(wfs(wf).adc_gains_dB(adc)/20) * chan_equal ...
        * 10.^(wfs(wf).system_dB/20)));
    else
      % A number is provided for each receiver path for system_dB
      mult_factor{img}(:,:,wf_adc) = single(wfs(wf).quantization_to_V(adc) ...
        / (10.^(wfs(wf).adc_gains_dB(adc)/20) * chan_equal ...
        * 10.^(wfs(wf).system_dB(param.radar.wfs(wf).rx_paths(adc))/20)));
    end
    
    if 1 %%% #######################################################################################
      raw_data{img}(:,:,wf_adc) = 1/mult_factor{img}(:,:,wf_adc) * raw_data{img}(:,:,wf_adc);
    end
    
  end %% for wf_adc
end %% for img


param.sim.mult_factor = mult_factor;
figure(h_fig); % for phase centers
% suptitle('Numbers = rx path; [O, *] = phase center [w/o, w] attitude  ');
set(findobj(h_fig,'type','axes'),'FontWeight', 'Bold', 'FontSize',14);
set(h_fig, 'Position', get(0, 'Screensize'));

%% Output_Files

false_save_en = 0; % if 1; does not save the files

% RECORDS
param.sim.out_fn_dir_records = fullfile(gRadar.support_path,'records', ...
  sprintf('%s',param.sim.radar_name), ...
  sprintf('%s',param.season_name) );
param.sim.out_fn_records = fullfile(param.sim.out_fn_dir_records, ...
  sprintf('records_%s.mat',param.day_seg) );
fprintf('(%s) Saving records \t%s\n', datestr(now), param.sim.out_fn_records);
if ~false_save_en; ct_save(param.sim.out_fn_records,'-struct', 'records'); end;

% FRAMES
param.sim.out_fn_dir_frames = fullfile(gRadar.support_path,'frames', ...
  sprintf('%s',param.sim.radar_name), ...
  sprintf('%s',param.season_name) );
param.sim.out_fn_frames = fullfile(param.sim.out_fn_dir_frames, ...
  sprintf('frames_%s.mat',param.day_seg) );
fprintf('(%s) Saving frames \t%s\n', datestr(now), param.sim.out_fn_frames);
if ~false_save_en; ct_save(param.sim.out_fn_frames, '-struct', 'frames'); end;

% RAW DATA
param.sim.out_fn_dir_raw_data = fullfile(gRadar.ct_tmp_path,'sim3D', ...
  sprintf('%s',param.sim.radar_name), ...
  sprintf('%s',param.season_name), ...
  sprintf('%s',param.day_seg(1:8)) );
param.sim.out_fn_raw_data = fullfile(param.sim.out_fn_dir_raw_data, ...
  sprintf('data_%s.mat', param.day_seg) );
fprintf('(%s) Saving raw_data \t%s\n', datestr(now), param.sim.out_fn_raw_data);
if ~false_save_en; ct_save(param.sim.out_fn_raw_data, 'raw_data');end;

% PARAM
param.sim.out_fn_dir_param = param.sim.out_fn_dir_raw_data;
param.sim.out_fn_param = fullfile(param.sim.out_fn_dir_param, ...
  sprintf('param.mat') );
param.fn = param.sim.out_fn_param;
fprintf('(%s) Saving param \t%s\n', datestr(now), param.sim.out_fn_param);
if ~false_save_en; ct_save(param.fn, 'param');end;

if false_save_en
  fprintf('(%s) ===== FALSE_SAVE_EN: none of the above files are written =====\n', datestr(now));
else
  fprintf('(%s) Files SAVED\n', datestr(now));
  fprintf('=====================================================================\n');
end

%% FIGURES (work best for single target case)

if strcmpi(param.target.type, 'point')
  point_target = 1;
else
  point_target = 0;
end

for compressing_this=1 %continue;
  
  fprintf('Max Power in each image [wf-adc](size)@(max_idx_row,max_idx_col))::dB (mul_factor_dB, adc_dB, system_dB)]: \n'); % per-image stats
  for img = 1:length(param.sim.imgs)
    for wf_adc = 1:size(param.sim.imgs{img},1)
      wf = param.sim.imgs{img}(wf_adc,1);
      adc = param.sim.imgs{img}(wf_adc,2);
      
      data = raw_data{img}(:,:,wf_adc);
      
      t_ref = param.radar.wfs(wf).time_raw;
      tmp = 20*log10(abs(data)) ;
      
      [tt_Nt, tt_Nx] = size(tmp);
      [tt_max, tt_idx] = max(tmp(:));
      [tt_row, tt_col] = ind2sub(size(tmp), tt_idx);
      if length(wfs(wf).system_dB) == 1
        tt_sys_dB = wfs(wf).system_dB;
      else
        tt_sys_dB = wfs(wf).system_dB(param.radar.wfs(wf).rx_paths(adc));
      end
      fprintf('\n[%d-%d](%d,%d)@(%d,%d)::%1.2f dB (mul:%1.2f,adc:%1.2fdB,sys:%1.2fdB)', ...
        wf, adc, tt_Nt, tt_Nx, tt_row, tt_col, tt_max, ...
        20*log10(mult_factor{img}(:,:,wf_adc)), wfs(wf).adc_gains_dB(adc),...
        tt_sys_dB ); % per-image stats
      
      %% figures
      
      if 0 || Nimgs > 9
        continue;
      end
      
      x = 1:Nx;
      if point_target
        [td_closest, idx_closest] = min(twtt);
      end
      
      call_sign = sprintf('Simulated Data wfs_%02d_adc_%02d',wf,adc);
      fig_title = sprintf('%s_%s',mfilename, call_sign);
      fig_h = figure('Name',fig_title);
      
      subplot(121);
      %imagesc(x, t_ref/1e-6, tmp-max(tmp(:)) ,[-30,0] );
      imagesc(x, t_ref/1e-6, tmp);
      cb = colorbar; cb.Label.String = 'Relative Power, dB';
      grid on; hold on; axis tight;
      if point_target
        plot(x, twtt/1e-6, '--', 'Color', 'g');
        plot(idx_closest, twtt(idx_closest)/1e-6, 's','LineWidth',2);
      end
      xlabel('Along-track position, rlines'); ylabel('Fast-time, us');
      title('PhaseHistory Magnitude');
      
      subplot(122);
      imagesc(x, t_ref/1e-6, angle(data) );
      cb = colorbar; cb.Label.String = 'Radians';
      grid on; hold on; axis tight;
      if point_target
        plot(x, twtt/1e-6, '--', 'Color', 'g');
        plot(idx_closest, twtt(idx_closest)/1e-6, 's','LineWidth',2);
      end
      xlabel('Along-track position, rlines'); ylabel('Fast-time, us');
      title('PhaseHistory Phase');
      
      try
        sgtitle(fig_title,'FontWeight','Bold','FontSize',14,'Interpreter','None');
      end
      set(findobj(fig_h,'type','axes'),'FontWeight', 'Bold', 'FontSize',14);
      set(fig_h, 'Position', get(0, 'Screensize'));
      %   print(gcf, '-dpng', fig_title, '-r300');
      
      if param.sim.debug_plots_en
        % Frequency check
        fig_title = sprintf('Frequency check wfs_%02d_adc_%02d',wf,adc);
        fig_h = figure('Name', fig_title);
        hold on;
        
        f0    = param.radar.wfs(wf).f0;
        f1    = param.radar.wfs(wf).f1;
        fs    = param.radar.fs;
        IF_nz = [0:3];
        IF_cutoffs = [0 : 0.5 : 2]' * fs ; % Hz
        
        for idx = 1:length(IF_cutoffs)-1
          rectangle('Position',[IF_cutoffs(idx)/1e6 0 fs/2e6 1.5],'EdgeColor','k');
        end
        rectangle('Position',[f0/1e6 0 (f1-f0)/1e6 1],'EdgeColor','r','FaceColor', 'g');
        xticks( round(sort( [IF_cutoffs; f0; f1] ./1e6 ),2) );
        yticks([1 1.5]);
        grid on; axis tight;
        ylim([0 2]);
        
        xlabel('Frequency, MHz'); ylabel('Frequency bands');
        title('4 Nyquist zones and passband');
      end
      
    end %% for wf_adc
    fprintf('\n'); % per-image stats
  end %% for img
  
end

return;
%% sanity checks

leg_str = [];
leg_str_diff = [];
leg_str_ref = [];
marker_char = [ '.' 'o' 's' '>' '<' '+' 'd' ];

pc = [];
pc_rel_ref = [];
pc_rel_ref_dist = [];
pc_to_target = [];
pc_to_target_dist = [];
pc_to_target_phase_delay = [];
diff_phase = [];

param.gps_source = records.gps_source;

for img = 1:length(param.sim.imgs)
  for wf_adc = 1:size(param.sim.imgs{img},1)
    wf = param.sim.imgs{img}(wf_adc,1);
    adc = param.sim.imgs{img}(wf_adc,2);
    
    leg_str = [leg_str {sprintf('[%d-%d]', wf, adc)} ];
    % param.radar.wfs(wf).rx_paths(adc)
    pc{img}(:,wf_adc) = [+1; -1; -1] .* lever_arm(param, [], param.radar.wfs(wf).rx_paths(adc) );
    
  end
end

for img = 1:length(param.sim.imgs)
  for wf_adc = 1:size(param.sim.imgs{img},1)
    wf = param.sim.imgs{img}(wf_adc,1);
    adc = param.sim.imgs{img}(wf_adc,2);
    
    if 0
      ref_img = 2;
      ref_wf_adc = 1; %%%%%%%%%%#########<<<<<<< just an idx,  not adc
    else
      ref_img = 1;
      ref_wf_adc = 2; %%%%%%%%%%#########<<<<<<< just an idx,  not adc
    end
    
    lambda = c / fc;
    lambda = c ./ [f0; fc; f1]; % to check phases at begining, center and end of phase ramp
    
    target = [ 0; 0; param.target.offsets.z ]; %%%%%%%%%%
    
    pc_rel_ref{img}(:,wf_adc) = pc{img}(:,wf_adc) - pc{ref_img}(:,ref_wf_adc) ;
    pc_rel_ref_dist{img}(wf_adc) = vecnorm(pc_rel_ref{img}(:,wf_adc));
    
    pc_to_target{img}(:,wf_adc) = target - pc_rel_ref{img}(:,wf_adc);
    pc_to_target_dist{img}(wf_adc) = vecnorm(pc_to_target{img}(:,wf_adc));
    
    % each col corresponds to a wf_adc for each row of [fo fc f1]
    pc_to_target_phase_delay{img}(:,wf_adc) = 2*2*pi./lambda * pc_to_target_dist{img}(wf_adc);
    
  end
  
  % each col corresponds to a wf_adc for each row of [fo fc f1]
  pc_diff_phase{img} = pc_to_target_phase_delay{img} - pc_to_target_phase_delay{img}(:,ref_wf_adc);
  
  t_ref = param.radar.wfs(wf).time_raw; % #####
  expected_beat_freq = (f1-f0)/Tpd * 2*(diff(pc_to_target_dist{img})) / c;
  dt = 1/param.radar.fs;
  t_dt = ( (t_ref+Tpd) - rem(t_ref,dt) ) ;
  t_dt = t_ref+Tpd;
  expected_phase_delay =  2*pi * expected_beat_freq .* t_dt ;
  
  %figure;
  % plot( t_ref, 2*2*pi * (f1-f0)/c * (pc_to_target_dist{1}(1)-pc_to_target_dist{1}(2)) * t_ref/Tpd)
end

for img = 1:length(param.sim.imgs)
  ttt = raw_data{img}(:,idx_closest,:);
  ttt = squeeze(ttt);
  figure;
  
  for wf_adc = 1:size(param.sim.imgs{img},1)
    wf = param.sim.imgs{img}(wf_adc,1);
    adc = param.sim.imgs{img}(wf_adc,2);
    
    t_ref = param.radar.wfs(wf).time_raw; % #####
    
    subplot(221);
    hold on;
    plot(t_ref/1e-6, 20*log10( abs(ttt(:,wf_adc)) ), [marker_char(wf_adc) '-'] );
    legend(leg_str{1:wf_adc});
    
    subplot(222);
    hold on;
    plot(t_ref/1e-6, angle( ttt(:,wf_adc) ), [marker_char(wf_adc) '-'] );
    legend(leg_str{1:wf_adc});
    
    subplot(223);
    if wf_adc>1
      hold on;
      plot(t_ref/1e-6, unwrap( angle(ttt(:,wf_adc)) - angle(ttt(:,wf_adc-1)) ), [marker_char(wf_adc) '-'] );
      leg_str_diff = [ leg_str_diff {sprintf('%s - %s',leg_str{wf_adc}, leg_str{wf_adc-1}) } ];
      legend( leg_str_diff );
    end
    
    
    subplot(224);
    if ref_wf_adc <= size(param.sim.imgs{img},1)
      hold on;
      plot(t_ref/1e-6, unwrap( angle(ttt(:,wf_adc)) - angle(ttt(:,ref_wf_adc)) ), [marker_char(wf_adc) '-'] );
      leg_str_ref = [ leg_str_ref { sprintf('%s - %s', leg_str{wf_adc}, leg_str{ref_wf_adc}) } ];
      legend( leg_str_ref );
    end
    
  end
end

subplot(223)
plot(t_ref/1e-6, expected_phase_delay, 's', 'Color', 'g');
leg_str_diff = [ leg_str_diff {'expected'} ];
legend( leg_str_diff );

subplot(224)
% plot(t_ref/1e-6, expected_phase_delay, 's', 'Color', 'g');
errorbar(t_ref/1e-6, expected_phase_delay(:,1),   2*pi * min(expected_beat_freq) * dt *ones(size(t_ref)) , 's', 'Color', 'g');
leg_str_ref = [ leg_str_ref {'expected'} ];
legend( leg_str_ref );

