
param = [];

%% Petermann Line 1 2014
% if ispc
%   fn = fullfile('X:\ct_data\rds\2014_Greenland_P3\CSARP_multipass\',sprintf('Petermann_line1_2014.mat'));
% else
%   fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_multipass/',sprintf('Petermann_line1_2014'));
% end

%% Petermann Line 1 2018
% if ispc
%   fn = fullfile('X:\ct_data\rds\2018_Greenland_P3\CSARP_multipass\',sprintf('Petermann_line1_2018.mat'));
% else
%   fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2018_Greenland_P3/CSARP_multipass/',sprintf('Petermann_line1_2018'));
% end

%% Petermann Line 1 2014, 2018
if ispc
  param.multipass.fn = fullfile('X:\ct_data\rds\2014_Greenland_P3\CSARP_multipass\',sprintf('Petermann_line1_2014_2018.mat'));
else
  param.multipass.fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_multipass/',sprintf('Petermann_line1_2014_2018'));
end

param.multipass.rbins = [];

param.multipass.baseline_master_idx = 1;
param.multipass.master_idx = 1;

param.multipass.pass_en_mask = [];
param.multipass.output_fn_midfix = [];
param.multipass.coregistration_time_shift = [0 -2];

%% Setup
% =========================================================================

proj = geotiffinfo(ct_filename_gis([],fullfile('greenland','Landsat-7','Greenland_natural_90m.tif')));

physical_constants;

%% Input check
% -----------------------

fn = param.multipass.fn;

if ~isfield(param.multipass,'baseline_master_idx') || isempty(param.multipass.baseline_master_idx)
  param.multipass.baseline_master_idx = 1;
end
baseline_master_idx = param.multipass.baseline_master_idx;

if ~isfield(param.multipass,'master_idx') || isempty(param.multipass.master_idx)
  param.multipass.master_idx = 1;
end
master_idx = param.multipass.master_idx;

if ~isfield(param.multipass,'output_fn_midfix') || isempty(param.multipass.output_fn_midfix)
  param.multipass.output_fn_midfix = '';
end
output_fn_midfix = param.multipass.output_fn_midfix;

if ~isfield(param.multipass,'units') || isempty(param.multipass.units)
  param.multipass.units = 'meters'; % 'meters' or 'bins'
end
master_idx = param.multipass.master_idx;

% Load multipass.combine_passes file
load(fn);

if ~isfield(param.multipass,'pass_en_mask') || isempty(param.multipass.pass_en_mask)
  param.multipass.pass_en_mask = [];
end
% Assume any undefined passes are enabled
param.multipass.pass_en_mask(end+1:length(pass)) = true;
pass_en_mask = param.multipass.pass_en_mask;

pass_en_idxs = find(pass_en_mask);

if ~isfield(param.multipass,'rbins') || isempty(param.multipass.rbins)
  param.multipass.rbins = [];
end
if isempty(param.multipass.rbins)
  param.multipass.rbins = 1:size(pass(baseline_master_idx).data,1);
end
rbins = param.multipass.rbins;
% HACK THAT NEEDS TO BE FIXED FOR MAKING ALL IMAGES SAME SIZE:
% for pass_idx = 1:length(pass)
%   if pass_en_mask(pass_idx)
%     rbins = intersect(rbins,1:size(pass(pass_idx).data,1));
%   end
% end

if ~isfield(param.multipass,'coregistration_time_shift') || isempty(param.multipass.coregistration_time_shift)
  param.multipass.coregistration_time_shift = zeros(1,length(pass));
end
coregistration_time_shift = param.multipass.coregistration_time_shift;

%% Plot Results
% =========================================================================
h_fig_map = figure(100); clf;
h_plot_map = [];
h_legend_map = {};
hold on;
axis('equal');
h_fig_elev = figure(101); clf;
h_plot_elev = [];
h_legend_elev = {};
hold on;
xlabel('Range line');
ylabel('WGS-84 elevation (m)');
grid on;

h_data_axes = [];
for pass_idx = 1:length(pass)
  if pass_en_mask(pass_idx)
    figure(pass_idx); clf;
    set(pass_idx,'WindowStyle','docked')
    imagesc(lp(pass(pass_idx).data(rbins,:)))
    colormap(1-gray(256));
    h_data_axes(end+1) = gca;
  end
  
  % Apply GPS time offset
  if 0
    pass_idx = 5;
    time_offset = -5;
    pass(pass_idx).ecef = interp1(pass(pass_idx).gps_time,pass(pass_idx).ecef.',pass(pass_idx).gps_time+time_offset,'linear','extrap').';
    pass(pass_idx).x = interp1(pass(pass_idx).gps_time,pass(pass_idx).x.',pass(pass_idx).gps_time+time_offset,'linear','extrap').';
  end
  
  pass(pass_idx).ecef = pass(pass_idx).origin;
  for rline = 1:size(pass(pass_idx).origin,2)
    pass(pass_idx).ecef(:,rline) = pass(pass_idx).ecef(:,rline) ...
      + [pass(pass_idx).x(:,rline) pass(pass_idx).y(:,rline) pass(pass_idx).z(:,rline)]*pass(pass_idx).pos(:,rline);
  end
  [pass(pass_idx).lat,pass(pass_idx).lon,pass(pass_idx).elev] = ecef2geodetic(referenceEllipsoid('wgs84'), ...
    pass(pass_idx).ecef(1,:), pass(pass_idx).ecef(2,:), pass(pass_idx).ecef(3,:));
  
  figure(h_fig_map);
  if 0
    h_plot = plot(pass(pass_idx).lon, pass(pass_idx).lat,'.');
    color = get(h_plot,'Color');
    h_text = text(pass(pass_idx).lon(1), pass(pass_idx).lat(1), sprintf('%d', pass_idx), 'Color', color);
  else
    [pass(pass_idx).proj_x,pass(pass_idx).proj_y] = projfwd(proj,pass(pass_idx).lat,pass(pass_idx).lon);
    h_plot_map(end+1) = plot(pass(pass_idx).proj_x/1e3, pass(pass_idx).proj_y/1e3,'.');
    h_legend_map{end+1} = sprintf('%d',pass_idx);
    color = get(h_plot_map(end),'Color');
    h_text = text(pass(pass_idx).proj_x(1)/1e3, pass(pass_idx).proj_y(1)/1e3, sprintf('%d', pass_idx), 'Color', color);
    xlabel('X (km)');
    ylabel('Y (km)');
    grid on;
  end
  
  figure(h_fig_elev);
  if baseline_master_idx == pass_idx
    h_plot_elev(end+1) = plot(pass(pass_idx).elev,'LineWidth',2,'UserData',pass_idx);
  else
    h_plot_elev(end+1) = plot(pass(pass_idx).elev,'UserData',pass_idx);
  end
  h_legend_elev{end+1} = sprintf('%d',pass_idx');
end
linkaxes(h_data_axes,'xy');
legend(h_plot_map,h_legend_map);

%% Co-Register Results
% =========================================================================


if 1
  % Option 1: Use a single pass as the reference.
  ref = pass(baseline_master_idx);
  
else
  % Option 2: Take a single pass as master and construct along track vectors relative
  % to it. Fit a polynomial to all the data using the along track vector
  % from step 1. Use this as the reference. This may be useful to do if no
  % single pass follows the center of the tube of passes.
  
  % With reference, run SAR coord system in a special mode where the mean is
  % not taken across Lsar, but instead each point is directly passed to the
  % output.
end

along_track = geodetic_to_along_track(ref.lat,ref.lon,ref.elev);
ref.surface_bin = interp1(ref.wfs(ref.wf).time, 1:length(ref.wfs(ref.wf).time), ref.surface);

if 0
  h_fig_ref_idx = figure(102); clf;
  hold on;
end

%% Pass
data = [];
for pass_idx = 1:length(pass)
  
  %% Pass: position in ref coordinate system
  pass(pass_idx).ref_idx = zeros(1,size(pass(pass_idx).origin,2));
  last_idx = 0;
  for rline = 1:size(pass(pass_idx).ecef,2)
    % offset: offset from position rline of pass pass_idx from every point
    %         on the ref line
    offset = bsxfun(@minus, pass(pass_idx).ecef(:,rline), ref.ecef);
    % dist: converts offset into along-track distance of the slave line
    dist = offset.'*pass(pass_idx).x(:,rline);
    % min_idx: finds the point on the reference line that is closest to
    %          this point
    [min_dist,min_idx] = min(abs(dist));
    %       if min_idx == last_idx
    %         keyboard
    %       end
    last_idx = min_idx;
    pass(pass_idx).ref_idx(rline) = min_idx;
    
    % x_offset: along-track offset on master line of the closest point
    x_offset = offset(:,min_idx).'*ref.x(:,min_idx);
    % Compute FCS of slave point in master line coordinate system
    %   FCS: flight (aka SAR) coordinate system
    pass(pass_idx).along_track(rline) = along_track(min_idx) + x_offset;
    pass(pass_idx).ref_y(rline) = offset(:,min_idx).'*ref.y(:,min_idx);
    pass(pass_idx).ref_z(rline) = offset(:,min_idx).'*ref.z(:,min_idx);
    
    if 0
      % Compute the location of all pixels from this range line in ECEF
      pass(pass_idx).wfs(pass(pass_idx).wf).time;
      time = pass(pass_idx).time(2)-pass(pass_idx).time(1);
      
      range = time * c/2;
      range(time>pass(pass_idx).surface(rline)) = pass(pass_idx).surface(rline)*c/2 ...
        + (time(time>pass(pass_idx).surface(rline)) - pass(pass_idx).surface(rline))*c/2/sqrt(er_ice);
      
      pixels = pass(pass_idx).ecef(:,rline) + pass(pass_idx).z(:,rline)*range;
      
      % Compute the location of all pixels from this range line in the master
      % line coordinate system FCS.
      Tfcs = [ref.x, ref.y, ref.z];
      pass(pass_idx).pixels(:,:,rline) = (pixels - ref.origin) / Tfcs;
    end
  end
  pass(pass_idx).along_track_slave = geodetic_to_along_track(pass(pass_idx).lat,pass(pass_idx).lon,pass(pass_idx).elev);
  
  if 0
    % Debug plot showing indexes for alignment of passes
    figure(h_fig_ref_idx);
    plot(pass(pass_idx).ref_idx)
    drawnow;
  end
  
  %% Pass: Resample in along-track
  % Resample images and position vectors onto a common along-track axes
  if 0
    % 1. Oversample slave data by 10x in along track
    Mx = 10;
    Nx = size(pass(pass_idx).data,2);
    data_oversample = interpft(pass(pass_idx).data.',Mx*Nx);
    % 2. Interpolate to find the oversampled slave axes
    along_track_oversample = interp1(0:Nx-1, ...
      pass(pass_idx).along_track, (0:Nx*Mx-1)/Mx,'linear','extrap');
    % 3. Interpolate oversampled slave data onto master along track axes
    pass(pass_idx).ref_data = interp1(along_track_oversample, ...
      data_oversample, along_track,'linear','extrap').';
  else
    pass(pass_idx).ref_data = interp1(pass(pass_idx).along_track, ...
      pass(pass_idx).data.', along_track,'linear','extrap').';
  end
  pass(pass_idx).ref_y = interp1(pass(pass_idx).along_track, ...
    pass(pass_idx).ref_y, along_track,'linear','extrap').';
  pass(pass_idx).ref_z = interp1(pass(pass_idx).along_track, ...
    pass(pass_idx).ref_z, along_track,'linear','extrap').';
  
  %% Pass: Apply fixed coregistration time shift
  Nt = size(pass(pass_idx).ref_data,1);
  dt = pass(pass_idx).wfs(pass(pass_idx).wf).time(2)-pass(pass_idx).wfs(pass(pass_idx).wf).time(1);
  time = dt*(0:Nt-1).';
  df = 1/(dt*Nt);
  freq = df * ifftshift( -floor(Nt/2) : floor((Nt-1)/2) ).';
%   freq = freq - freq(1); % Remove center frequency offset
  dt = coregistration_time_shift(pass_idx) * (pass(pass_idx).wfs(pass(pass_idx).wf).time(2)-pass(pass_idx).wfs(pass(pass_idx).wf).time(1));
  if 0
    pass(pass_idx).ref_data = ifft(bsxfun(@times,fft(double(pass(pass_idx).ref_data)),exp(-1i*2*pi*freq*dt)));
  else
    pass(pass_idx).ref_data = interp1(pass(pass_idx).wfs(pass(pass_idx).wf).time, pass(pass_idx).ref_data, pass(pass_idx).wfs(pass(pass_idx).wf).time+dt, 'linear');
    pass(pass_idx).ref_data = interp_finite(pass(pass_idx).ref_data);
  end
  
  %% Pass: Motion/slope compensation
  insar_mode = 2; % HACK!!!
  if insar_mode == 1
    % Motion compensation of FCS z-motion
    for rline = 1:size(pass(pass_idx).ref_data,2)
      % Convert z-offset into time-offset assuming nadir DOA
      dt = pass(pass_idx).ref_z(rline)/(c/2);
      pass(pass_idx).ref_data(:,rline) = ifft(fft(pass(pass_idx).ref_data(:,rline)) ...
        .*exp(1i*2*pi*pass(pass_idx).wfs(pass(pass_idx).wf).freq*dt) );
    end
  elseif insar_mode == 2 || insar_mode == 4
    % Co-register images using GPS and nadir squint angle assumption
    %
    if 0
      % Motion compensation of FCS z-motion without center frequency so there
      % is no phase shift.
      %     freq = pass(pass_idx).wfs(pass(pass_idx).wf).freq;
      %     freq = freq - freq(1); % Remove center frequency offset
      for rline = 1:size(pass(pass_idx).ref_data,2)
        % Convert z-offset into time-offset assuming nadir DOA
        dt = pass(pass_idx).ref_z(rline)/(c/2);
        pass(pass_idx).ref_data(:,rline) = ifft(fft(pass(pass_idx).ref_data(:,rline)) ...
          .*exp(1i*2*pi*freq*dt) );
      end
    else
      % Motion compensation of FCS z-motion using linear interpolation
      for rline = 1:size(pass(pass_idx).ref_data,2)
        dt = pass(pass_idx).ref_z(rline)/(c/2);
        pass(pass_idx).ref_data(:,rline) = interp1(pass(pass_idx).wfs(pass(pass_idx).wf).time, pass(pass_idx).ref_data(:,rline), pass(pass_idx).wfs(pass(pass_idx).wf).time+dt, 'linear');
        pass(pass_idx).ref_data(:,rline) = interp_finite(pass(pass_idx).ref_data(:,rline));
      end
    end
  elseif insar_mode == 3
    % Motion compensation of FCS z-motion and slope compensation
    
    % True time delay shift for z-offset
    for rline = 1:size(pass(pass_idx).ref_data,2)
      % Convert z-offset into time-offset assuming nadir DOA
      dt = pass(pass_idx).ref_z(rline)/(c/2);
      pass(pass_idx).ref_data(:,rline) = ifft(fft(pass(pass_idx).ref_data(:,rline)) ...
        .*exp(1i*2*pi*pass(pass_idx).wfs(pass(pass_idx).wf).freq*dt) );
    end
    
    % Phase only correction for slope
    if 1
      % Using file generated from this dataset
      [fn_dir,fn_name] = fileparts(fn);
      fn_slope = fullfile(fn_dir,[fn_name '_slope.mat']);
      load(fn_slope,'slope','GPS_time','Latitude','Longitude','Elevation','Time','Surface');
      slope = interp1(GPS_time,slope.',pass(baseline_master_idx).gps_time).';
      slope = interp_finite(slope.').';
      slope = interp1(Time,slope,pass(pass_idx).wfs(pass(pass_idx).wf).time);
      slope = interp_finite(slope);
      
      pass(pass_idx).ref_data = pass(pass_idx).ref_data .* exp(-1i*4*pi*pass(pass_idx).wfs(pass(pass_idx).wf).freq(1)/c *bsxfun(@times,sin(slope),pass(pass_idx).ref_y(:).'));
      
    elseif 1
      % Using file generated from another dataset
      fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_20140429_01_067_wf2_slope.mat';
      
      % TBD
      
    end
    
  elseif 0
    % Co-register images using cross-correlation
    keyboard
  end
  
  %% Pass: Match time axis to baseline_master_idx
  if 0
    pass(pass_idx).ref_data = interp1(pass(pass_idx).wfs(pass(pass_idx).wf).time, pass(pass_idx).ref_data, pass(baseline_master_idx).wfs(pass(baseline_master_idx).wf).time, 'linear', 0);
  else
    Mt = 4;
    Nt = length(pass(pass_idx).wfs(pass(pass_idx).wf).time);
    dt = pass(pass_idx).wfs(pass(pass_idx).wf).time(2)-pass(pass_idx).wfs(pass(pass_idx).wf).time(1);
    if 0
      pass(pass_idx).ref_data = interpft(pass(pass_idx).ref_data,Mt*Nt);
      time_Mt = pass(pass_idx).wfs(pass(pass_idx).wf).time(1) + dt/Mt*(0:Mt*Nt-1);
      pass(pass_idx).ref_data = interp1(time_Mt, pass(pass_idx).ref_data, pass(baseline_master_idx).wfs(pass(baseline_master_idx).wf).time, 'linear', 0);
    else
      pass(pass_idx).ref_data = interp1(pass(pass_idx).wfs(pass(pass_idx).wf).time, pass(pass_idx).ref_data, pass(baseline_master_idx).wfs(pass(baseline_master_idx).wf).time, 'linear', 0);
    end
  end
  
  if 0
    % Normalize surface phase
    Nt = size(pass(pass_idx).ref_data,1);
    Nx = size(pass(pass_idx).ref_data,2);
    H = pass(baseline_master_idx).ref_data(round(ref.surface_bin)+(0:Nx-1)*Nt) .* conj(pass(pass_idx).ref_data(round(ref.surface_bin)+(0:Nx-1)*Nt));
    H = exp(1i*angle(H));
    pass(pass_idx).ref_data = bsxfun(@times,pass(pass_idx).ref_data,H);
  end
  
  % Concatenate data into a single matrix
  data = cat(3,data,single(pass(pass_idx).ref_data));
end

%% Apply equalization
% -----------------------
% if insar_mode == 2 || insar_mode == 3 || insar_mode == 4
%   equalization = reshape(equalization,[1 1 numel(equalization)]);
%   data(:,:,pass_en_idxs) = bsxfun(@times,data(:,:,pass_en_idxs),1./equalization(:,:,pass_en_idxs));
% end

if 0
  %% Coregister: Data Dependent method to estimate System Time Delay
  % Apply fixed coregistration time shift
  for pass_out_idx = 2%1:length(pass_en_idxs)
    pass_idx = pass_en_idxs(pass_out_idx);
    freq = pass(pass_idx).wfs(pass(pass_idx).wf).freq;
    freq = freq - freq(1); % Remove center frequency offset
    coregistration_time_shifts = -2:0.05:2;
    coregistration_time_shifts = -1.6:0.01:-1.3;
    %     coregistration_time_shifts = -0.2:0.01:0.2;
    coherence_sum = [];
    for coregistration_time_shift_idx = 1:length(coregistration_time_shifts)
      coregistration_time_shift = coregistration_time_shifts(coregistration_time_shift_idx);
      dt = coregistration_time_shift * (pass(pass_idx).wfs(pass(pass_idx).wf).time(2)-pass(pass_idx).wfs(pass(pass_idx).wf).time(1));
      adjusted = ifft(bsxfun(@times,fft(data(:,:,pass_idx)),exp(-1i*2*pi*freq*dt)));
      coherence = fir_dec(adjusted(rbins,:) .* conj(data(rbins,:,master_idx)) ./ abs(adjusted(rbins,:) .* data(rbins,:,master_idx)),ones(1,7)/7,1);
      coherence = fir_dec(coherence.',ones(1,3)/3,1).';
      coherence = abs(coherence);
      %     coherence_sum(coregistration_time_shift_idx) = sum(coherence(coherence>0.5));
      coherence_sum(coregistration_time_shift_idx) = sum(coherence(coherence>0));
      %    TriangleRayIntersection fprintf('%g %.2f\n', coregistration_time_shift, coherence_sum(coregistration_time_shift_idx));
      %     imagesc(coherence); colormap(1-gray(256));
      %     pause
    end
  end
  figure(1000); clf;
  plot(coregistration_time_shifts,coherence_sum)
  [~,coregistration_time_shift_idx] = max(coherence_sum);
  coregistration_time_shift = coregistration_time_shifts(coregistration_time_shift_idx)
  return
end

%% Plot interferograms
h_data_axes = [];
new_equalization = [];
for pass_out_idx = 1:length(pass_en_idxs)
  pass_idx = pass_en_idxs(pass_out_idx);
  
  figure(pass_idx); clf;
  set(pass_idx,'WindowStyle','docked')
  if 1
    imagesc(lp(data(rbins,:,pass_idx)))
    colormap(1-gray(256));
    ylabel('Range bin');
    xlabel('Range line');
    title(sprintf('%s_%03d %d',pass(pass_idx).param_sar.day_seg,pass(pass_idx).frms(1),pass(pass_idx).direction),'interpreter','none')
    %caxis([-90 8]);
  else
    % Form interferogram (couple options)
    complex_data = fir_dec(data(rbins,:,pass_idx) .* conj(data(rbins,:,master_idx)),ones(1,11)/11,1);
    if ~exist('equalization_rlines','var') || isempty(equalization_rlines)
      new_equalization(pass_idx) = mean(complex_data(:)); % equalization only valid when motion compensation with phase is used
    else
      new_equalization(pass_idx) = mean(mean(complex_data(:,equalization_rlines))); % equalization only valid when motion compensation with phase is used
    end
    if insar_mode == 1
      complex_data = complex_data ./ new_equalization(pass_idx);
    end
    % Plot interferogram
    if insar_mode == 4
      imagesc(lp(data(rbins,:,pass_idx)));
      colormap(1-gray(256));
      h_colorbar = colorbar;
      set(get(h_colorbar,'ylabel'),'string','Relative power (dB)');
    else
      coherence = abs(fir_dec(data(rbins,:,pass_idx) .* conj(data(rbins,:,master_idx)) ./ abs(data(rbins,:,pass_idx) .* data(rbins,:,master_idx)),ones(1,11)/11,1)) ...
        .* exp(1i*angle(complex_data));
      imagesc(hsv_plot_coherence(coherence,[0 1]));
      colormap(hsv(256))
      h_colorbar = colorbar;
      caxis([-pi pi])
      set(get(h_colorbar,'ylabel'),'string','Angle (radians)');
    end
    ylabel('Range bin');
    xlabel('Range line');
    title(sprintf('%s_%03d %d',pass(pass_idx).param_sar.day_seg,pass(pass_idx).param_sar.load.frm,pass(pass_idx).direction),'interpreter','none')
    
    if 0
      imagesc(abs(coherence));
      colormap(1-gray(256))
      h_colorbar = colorbar;
      caxis([0 1])
      set(get(h_colorbar,'ylabel'),'string','Coherence');
      ylabel('Range bin');
      xlabel('Range line');
      title(sprintf('%s_%03d %d',pass(pass_idx).param_sar.day_seg,pass(pass_idx).param_sar.load.frm,pass(pass_idx).direction),'interpreter','none')
    end
    
  end
  h_data_axes(end+1) = gca;
  
end
% fprintf('=============================================\n');
% fprintf('New equalization\n');
% fprintf('%.1f ', lp(new_equalization)-mean(lp(new_equalization(pass_en_idxs))));
% fprintf('\n');
% fprintf('%.1f ', angle(new_equalization)*180/pi)
% fprintf('\n');
% fprintf('=============================================\n');
linkaxes(h_data_axes,'xy');
if insar_mode == 1 || insar_mode == 4
  return
end
if insar_mode == 2
  [fn_dir,fn_name] = fileparts(fn);
  fn_multipass = fullfile(fn_dir,[fn_name '_multipass.mat']);
  param_sar = pass(master_idx).param_sar;
  param_records = pass(master_idx).param_records;
  save(fn_multipass,'-v7.3','data','ref','param_sar','param_records');
  return
end
