
%% Setup
% =========================================================================
fprintf('=====================================================================\n');
fprintf('%s [Mode %d]: %s  (%s)\n', mfilename, param.multipass.comp_mode, ...
  param.multipass.pass_name, datestr(now));
fprintf('=====================================================================\n');
physical_constants;
standard_projections;

%% Load multipass.combine_passes file
fn = param.multipass.fn;

load(fn);  

%% Input check
% =========================================================================
% Confirm either SAR or echogram data
if ~isfield(param.multipass,'echo_sar') 
  if ~isfield(pass(1),'echo_sar')
      param.multipass.echo_sar = 'echo'; 
  else
    param.multipass.echo_sar = pass(1).echo_sar;
  end
end

% Make equalization
if ~isfield(param.multipass, 'equalization')
  param.multipass.equalization = exp(1i*(zeros(1,length(pass))/20)/180*pi);
end
equalization = param.multipass.equalization;
% All images are registered to the pass indicated by the baseline_master_idx
% baseline_master_idx does not need to be enabled
if ~isfield(param.multipass,'baseline_master_idx') || isempty(param.multipass.baseline_master_idx)
  param.multipass.baseline_master_idx = 1;
end
baseline_master_idx = param.multipass.baseline_master_idx;

if ~isfield(param.multipass,'coregistration_time_shift') || isempty(param.multipass.coregistration_time_shift)
  param.multipass.coregistration_time_shift = zeros(1,length(pass));
end
if length(param.multipass.coregistration_time_shift) < length(pass)
  param.multipass.coregistration_time_shift(length(pass)) = 0;
end
coregistration_time_shift = param.multipass.coregistration_time_shift;

if ~isfield(param.multipass,'debug_plots') || isempty(param.multipass.debug_plots)
  param.multipass.debug_plots = {'debug'};
end
enable_debug_plot = any(strcmp('debug',param.multipass.debug_plots));

if ~isfield(param.multipass,'layer') || isempty(param.multipass.layer)
  param.multipass.layer = struct();
end
if length(param.multipass.layer) < 2
  param.multipass.layer(2) = struct();
end
if ~isfield(param.multipass.layer(1),'name') || isempty(param.multipass.layer(1).name)
  param.multipass.layer(1).name = 'surface';
end
if ~isfield(param.multipass.layer(1),'source') || isempty(param.multipass.layer(1).source)
  param.multipass.layer(1).source = 'layerData';
end
if ~isfield(param.multipass.layer(2),'name') || isempty(param.multipass.layer(2).name)
  param.multipass.layer(2).name = 'bottom';
end
if ~isfield(param.multipass.layer(1),'source') || isempty(param.multipass.layer(2).source)
  param.multipass.layer(2).source = 'layerData';
end

% All comparisons are done relative to the pass indicated by the master_idx
% master_idx must be enabled
if ~isfield(param.multipass,'master_idx') || isempty(param.multipass.master_idx)
  param.multipass.master_idx = 1;
end
master_idx = param.multipass.master_idx;

% output_fn_midfix: string to add to output filenames (default is empty string)
if ~isfield(param.multipass,'output_fn_midfix') || isempty(param.multipass.output_fn_midfix)
  param.multipass.output_fn_midfix = '';
end
output_fn_midfix = param.multipass.output_fn_midfix;

% pass_en_mask: logical mask for each pass (default is all true)
if ~isfield(param.multipass,'pass_en_mask') || isempty(param.multipass.pass_en_mask)
  param.multipass.pass_en_mask = [];
end
% Assume any undefined passes are enabled
param.multipass.pass_en_mask(end+1:length(pass)) = true;
pass_en_mask = param.multipass.pass_en_mask;
if ~pass_en_mask(master_idx)
  error('The master index pass %d must be enabled in param.multipass.pass_en_mask.', master_idx);
end
pass_en_idxs = find(pass_en_mask);

if strcmpi(pass(baseline_master_idx).param_multipass.post.ops.location,'arctic')
  proj = arctic_proj;
elseif strcmpi(pass(baseline_master_idx).param_multipass.post.ops.location,'antarctic')
  proj = antarctic_proj;
else
  error('Unsupported location pass(%d).param_multipass.post.ops.location.', baseline_master_idx);
end

if ~isfield(param.multipass,'time_gate') || isempty(param.multipass.time_gate)
  param.multipass.time_gate = [-inf inf];
end

if ~isfield(param.multipass,'units') || isempty(param.multipass.units)
  param.multipass.units = 'meters'; % 'meters' or 'bins'
end

%% Convert FCS to ECEF and Geodetic
% =========================================================================
min_twtt = inf;
max_twtt = -inf;
max_rlines = -inf;
for pass_idx = 1:length(pass)
  if ~pass_en_mask(pass_idx) && pass_idx ~= baseline_master_idx
    pass(pass_idx).data = []; % Save some memory
    continue;
  end
  
  % Debug: GPS offset
  if 0
    pass_idx = 5;
    time_offset = -5;
    pass(pass_idx).ecef = interp1(pass(pass_idx).gps_time,pass(pass_idx).ecef.',pass(pass_idx).gps_time+time_offset,'linear','extrap').';
    pass(pass_idx).x = interp1(pass(pass_idx).gps_time,pass(pass_idx).x.',pass(pass_idx).gps_time+time_offset,'linear','extrap').';
  end
  
  % Convert SAR coordinate system (aka flight coordinate system, FCS) to
  % earth centered earth fixed (ECEF)
  pass(pass_idx).ecef = pass(pass_idx).origin;
  for rline = 1:size(pass(pass_idx).origin,2)
    pass(pass_idx).ecef(:,rline) = pass(pass_idx).ecef(:,rline) ...
      + [pass(pass_idx).x(:,rline) pass(pass_idx).y(:,rline) pass(pass_idx).z(:,rline)]*pass(pass_idx).pos(:,rline);
  end
  
  % Convert ECEF to Geodetic
  [pass(pass_idx).lat,pass(pass_idx).lon,pass(pass_idx).elev] = ecef2geodetic(referenceEllipsoid('wgs84'), ...
    pass(pass_idx).ecef(1,:), pass(pass_idx).ecef(2,:), pass(pass_idx).ecef(3,:));
  
  % Convert Geodetic to projected coordinates
  [pass(pass_idx).proj_x,pass(pass_idx).proj_y] = projfwd(proj,pass(pass_idx).lat,pass(pass_idx).lon);
  
  min_twtt = min(min_twtt,pass(pass_idx).wfs(pass(pass_idx).wf).time(1));
  max_twtt = max(max_twtt,pass(pass_idx).wfs(pass(pass_idx).wf).time(end));
  max_rlines = max(max_rlines,length(pass(pass_idx).gps_time));
end
param.multipass.time_gate(1) = max(param.multipass.time_gate(1),min_twtt);
param.multipass.time_gate(2) = min(param.multipass.time_gate(2),max_twtt);

%% Load layers
% =========================================================================
for pass_idx = 1:length(pass)
  if ~pass_en_mask(pass_idx)
    continue
  end
  if ~isfield(pass,'layers')
    % Load layers
    pass(pass_idx).layers = opsLoadLayers(pass(pass_idx).param_multipass,param.multipass.layer);
  end

  % Interpolate all layers onto a common reference (ref)
  for lay_idx = 1:length(pass(pass_idx).layers)
    pass(pass_idx).layers(lay_idx).twtt ...
      = interp_finite(interp1(pass(pass_idx).layers(lay_idx).gps_time, ...
      pass(pass_idx).layers(lay_idx).twtt, ...
      pass(pass_idx).gps_time,'linear'));
  end
end

%% Plot Raw Passes
% =========================================================================
if enable_debug_plot
  h_fig_map = figure(1001); clf(h_fig_map);
  h_plot_map = [];
  h_legend_map = {};
  h_axes_map = axes('parent',h_fig_map);
  hold(h_axes_map,'on');
  axis(h_axes_map, 'equal');
  xlabel(h_axes_map, 'X (km)');
  ylabel(h_axes_map, 'Y (km)');
  grid(h_axes_map, 'on');
  
  h_fig_elev = figure(1002); clf(h_fig_elev);
  h_plot_elev = [];
  h_legend_elev = {};
  h_axes_elev = axes('parent',h_fig_elev);
  hold(h_axes_elev,'on');
  xlabel(h_axes_elev, 'Range line');
  ylabel(h_axes_elev, 'WGS-84 elevation (m)');
  grid(h_axes_elev, 'on');
  
  clear h_axes_echo;
  for pass_idx = 1:length(pass)
    if ~pass_en_mask(pass_idx)
      continue
    end
    
    %% Plot: 1. plot echogram
    h_fig_echo = figure(pass_idx); clf(h_fig_echo);
    set(h_fig_echo,'WindowStyle','docked')
    set(h_fig_echo,'NumberTitle','off')
    set(h_fig_echo,'Name',num2str(pass_idx))
    h_axes_echo(pass_idx) = axes('parent',h_fig_echo);
    imagesc([],pass(pass_idx).wfs(pass(pass_idx).wf).time*1e6,lp(pass(pass_idx).data),'parent', h_axes_echo(pass_idx));
    title_str = pass(pass_idx).param_multipass.day_seg;
    title_str = regexprep(title_str,'_','\\_');
    title(h_axes_echo(pass_idx),title_str);
    colormap(h_axes_echo(pass_idx), 1-gray(256));
    hold(h_axes_echo(pass_idx),'on');
    if 1
      xlabel(h_axes_echo(pass_idx), 'Range line');
      ylabel(h_axes_echo(pass_idx), 'Two way travel time (\mus)');
    else
      [h_axes_echo_background,hp1,hp2] = plotyy(0:Nt-1,0:Nt-1,pass(pass_idx).wfs(pass(pass_idx).wf).time*1e6,pass(pass_idx).wfs(pass(pass_idx).wf).time*1e6,'parent',h_fig_echo);
      ylim(h_axes_echo_background(1),[0 Nt-1]);
      ylim(h_axes_echo_background(2),pass(pass_idx).wfs(pass(pass_idx).wf).time([1 end])*1e6);
      xlabel(h_axes_echo_background(1), 'Range line');
      ylabel(h_axes_echo_background(1), 'Range bin');
      ylabel(h_axes_echo_background(2), 'Time (\mus)');
      set(hp1,'XData',NaN,'YData',NaN);
      set(hp2,'XData',NaN,'YData',NaN);
      set(h_axes_echo_background(2),'YDir','reverse')
    end
    
    %% Plot: 2. plot layers
    Nx = size(pass(pass_idx).data,2);
    for lay_idx = 1:length(pass(pass_idx).layers)
      plot(1:Nx,pass(pass_idx).layers(lay_idx).twtt*1e6,'parent',h_axes_echo(pass_idx));
    end
    
    %% Plot: 3. plot map
    h_plot_map(end+1) = plot(pass(pass_idx).proj_x/1e3, pass(pass_idx).proj_y/1e3,'.','parent',h_axes_map);
    h_legend_map{end+1} = sprintf('%d',pass_idx);
    color = get(h_plot_map(end),'Color');
    h_text = text(pass(pass_idx).proj_x(1)/1e3, pass(pass_idx).proj_y(1)/1e3, sprintf('%d', pass_idx), 'Color', color,'parent',h_axes_map);
    
    %% Plot: 4. plot elevation
    if baseline_master_idx == pass_idx
      h_plot_elev(end+1) = plot(pass(pass_idx).elev,'LineWidth',2,'UserData',pass_idx,'parent',h_axes_elev);
    else
      h_plot_elev(end+1) = plot(pass(pass_idx).elev,'UserData',pass_idx,'parent',h_axes_elev);
    end
    h_legend_elev{end+1} = sprintf('%d',pass_idx');
  end
  legend(h_plot_map,h_legend_map);
  legend(h_plot_elev,h_legend_elev);
  linkaxes(h_axes_echo);
  xlim(h_axes_echo(1),[1 max_rlines]);
  ylim(h_axes_echo(1),param.multipass.time_gate*1e6);
end

%% Define reference path
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

if 0
  h_fig_ref_idx = figure(102); clf;
  hold on;
end

%% Pass Processing
% =========================================================================
data = [];
surf_flatten_en = false;
if surf_flatten_en
  ref.surface_bin = interp1(ref.wfs(ref.wf).time, 1:length(ref.wfs(ref.wf).time), ref.surface);
end
for pass_idx = 1:length(pass)
  %% Pass: 1. Position in ref coordinate system
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
    
    if 0 %strcmp('sar',param.multipass.echo_sar)
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
  
  %% Pass: 2. Resample in along-track
  % Resample images and position vectors onto a common along-track axes
  if strcmp('sar',param.multipass.echo_sar)
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
  elseif strcmp('echo',param.multipass.echo_sar)
    pass(pass_idx).ref_data = interp1(pass(pass_idx).along_track, ...
      pass(pass_idx).data.', along_track,'linear').';
  end
  
  pass(pass_idx).ref_y = interp1(pass(pass_idx).along_track, ...
    pass(pass_idx).ref_y, along_track,'linear','extrap').';
  pass(pass_idx).ref_z = interp1(pass(pass_idx).along_track, ...
    pass(pass_idx).ref_z, along_track,'linear','extrap').';
  
  for lay_idx = 1:length(pass(pass_idx).layers)
    pass(pass_idx).layers(lay_idx).twtt_ref = interp1(pass(pass_idx).along_track, ...
      pass(pass_idx).layers(lay_idx).twtt, along_track,'linear').';
  end
  
  %% Pass: 3. Apply fixed coregistration time shift
  if strcmp('echo',param.multipass.echo_sar)
    Nt = size(pass(pass_idx).ref_data,1);
    dt = pass(pass_idx).wfs(pass(pass_idx).wf).time(2)-pass(pass_idx).wfs(pass(pass_idx).wf).time(1);
    time = dt*(0:Nt-1).';
    df = 1/(dt*Nt);
    freq = df * ifftshift( -floor(Nt/2) : floor((Nt-1)/2) ).';
  elseif strcmp('sar',param.multipass.echo_sar)
    freq = []; time = [];
    for wf_id = 1:length(pass(pass_idx).wfs)
      freq = [freq pass(pass_idx).wfs(wf_id).freq'];
      time = [time pass(pass_idx).wfs(wf_id).time'];
    end
    %Sort the time values to be monotonically increasing
    [time, timesortind] =sort(time);
    freq = freq(timesortind);
    %Get rid of repeating values
    [time, timeunqind] = unique(time);
    freq = freq(timeunqind);
    freq = interp1(time,freq,pass(pass_idx).time);
    if param.multipass.comp_mode == 2
      % Motion compensation of FCS z-motion without center frequency so there
      freq = freq - freq(1); % Remove center frequency offset
    end
    dt = coregistration_time_shift(pass_idx) * (pass(pass_idx).wfs(pass(pass_idx).wf).time(2)-pass(pass_idx).wfs(pass(pass_idx).wf).time(1));
  end
  
  if strcmp('sar',param.multipass.echo_sar)
%     pass(pass_idx).ref_data = ifft(bsxfun(@times,fft(double(pass(pass_idx).ref_data)),exp(-1i*2*pi*freq*dt)));
    pass(pass_idx).ref_data = ifft(bsxfun(@times,fft(pass(pass_idx).ref_data),exp(-1i*2*pi*freq*dt)));
  elseif strcmp('echo',param.multipass.echo_sar)
    pass(pass_idx).ref_data = interp1(pass(pass_idx).wfs(pass(pass_idx).wf).time, pass(pass_idx).ref_data, pass(pass_idx).wfs(pass(pass_idx).wf).time+dt, 'linear');
    pass(pass_idx).ref_data = interp_finite(pass(pass_idx).ref_data);
    for lay_idx = 1:length(pass(pass_idx).layers)
      pass(pass_idx).layers(lay_idx).twtt_ref = pass(pass_idx).layers(lay_idx).twtt_ref - dt;
    end
  end
  
  %% Pass: 4. Motion/slope compensation
  if strcmp('sar',param.multipass.echo_sar) || any(param.multipass.comp_mode==[1 3 4])
    % Motion compensation of FCS z-motion
    for rline = 1:size(pass(pass_idx).ref_data,2)
      % Convert z-offset into time-offset assuming nadir DOA
      dt = pass(pass_idx).ref_z(rline)/(c/2);
      pass(pass_idx).ref_data(:,rline) = ifft(fft(pass(pass_idx).ref_data(:,rline)) ...
        .*exp(1i*2*pi*freq*dt) );
    end
  elseif strcmp('echo',param.multipass.echo_sar) && param.multipass.comp_mode == 2
    % Motion compensation of FCS z-motion using linear interpolation
    for rline = 1:size(pass(pass_idx).ref_data,2)
      dt = pass(pass_idx).ref_z(rline)/(c/2);
      pass(pass_idx).ref_data(:,rline) = interp1(pass(pass_idx).wfs(pass(pass_idx).wf).time, pass(pass_idx).ref_data(:,rline), pass(pass_idx).wfs(pass(pass_idx).wf).time+dt, 'linear');
      pass(pass_idx).ref_data(:,rline) = interp_finite(pass(pass_idx).ref_data(:,rline));
      for lay_idx = 1:length(pass(pass_idx).layers)
        pass(pass_idx).layers(lay_idx).twtt_ref(rline) = pass(pass_idx).layers(lay_idx).twtt_ref(rline) - dt;
      end
    end
  end
  
  if param.multipass.comp_mode == 3    
    % Phase only correction for slope
    if 1
      % Using file generated from this dataset
      [fn_dir,fn_name] = fileparts(fn);
      fn_slope = fullfile(fn_dir,[fn_name '_slope.mat']);
      load(fn_slope,'slope','GPS_time','Latitude','Longitude','Elevation','Time','Surface');
      slope = interp1(GPS_time,slope.',pass(baseline_master_idx).gps_time).';
      slope = interp_finite(slope.',0).';
      slope = interp1(Time,slope,pass(pass_idx).wfs(pass(pass_idx).wf).time);
      slope = interp_finite(slope);
      
      pass(pass_idx).ref_data = pass(pass_idx).ref_data .* exp(-1i*4*pi*freq(1)/c *bsxfun(@times,sin(slope),pass(pass_idx).ref_y(:).'));
      
    elseif 0
      % Using file generated from another dataset
      fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_20140429_01_067_wf2_slope.mat';
      
      % TBD
      
    end
    
  elseif 0
    % Co-register images using cross-correlation
    keyboard
  end
  
  %% Pass: Match time axis to baseline_master_idx
  % =========================================================================
  if 0
    pass(pass_idx).ref_data = interp1(pass(pass_idx).wfs(pass(pass_idx).wf).time, pass(pass_idx).ref_data, pass(baseline_master_idx).wfs(pass(baseline_master_idx).wf).time, 'linear', 0);
  elseif 1
    Mt = 4;
    Nt = length(pass(pass_idx).wfs(pass(pass_idx).wf).time);
    dt = pass(pass_idx).wfs(pass(pass_idx).wf).time(2)-pass(pass_idx).wfs(pass(pass_idx).wf).time(1);
    if strcmp('sar',param.multipass.echo_sar)
      pass(pass_idx).ref_data = interpft(pass(pass_idx).ref_data,Mt*Nt);
      time_Mt = pass(pass_idx).wfs(pass(pass_idx).wf).time(1) + dt/Mt*(0:Mt*Nt-1);
      pass(pass_idx).ref_data = interp1(time_Mt, pass(pass_idx).ref_data, pass(baseline_master_idx).wfs(pass(baseline_master_idx).wf).time, 'linear', 0);
    elseif strcmp('echo',param.multipass.echo_sar)
      pass(pass_idx).ref_data = interp1(pass(pass_idx).wfs(pass(pass_idx).wf).time, pass(pass_idx).ref_data, pass(baseline_master_idx).wfs(pass(baseline_master_idx).wf).time, 'linear', 0);
    end
  end
  
  if surf_flatten_en
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
if param.multipass.comp_mode ~= 1 && strcmp('sar',param.multipass.echo_sar)
  equalization = reshape(equalization,[1 1 numel(equalization)]);
  data(:,:,pass_en_idxs) = bsxfun(@times,data(:,:,pass_en_idxs),1./equalization(:,:,pass_en_idxs));
end

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
  figure(1002); clf;
  plot(coregistration_time_shifts,coherence_sum)
  [~,coregistration_time_shift_idx] = max(coherence_sum);
  coregistration_time_shift = coregistration_time_shifts(coregistration_time_shift_idx)
  return
end

%% Plot interferograms
h_data_axes = [];
new_equalization = [];
rbins = 1:size(data,1);
for pass_out_idx = 1:length(pass_en_idxs)
  pass_idx = pass_en_idxs(pass_out_idx);
  
  figure(pass_idx); clf;
  set(pass_idx,'WindowStyle','docked')
  if 1
    if strcmp(param.multipass.units,'meters')
      img = lp(data(rbins,:,pass_idx));
      
      elevation = pass(master_idx).elev;
      time = pass(master_idx).wfs(pass(pass_idx).wf).time(rbins);
      surface = pass(pass_idx).layers(1).twtt_ref(:).';
      %surface = pass(master_idx).surface;
      elev_max = max(elevation - time(1)*c/2);
      elev_min = min(elevation - surface*c/2 - (time(end)-surface)*c/(sqrt(er_ice)*2));
      dt = time(2) - time(1);
      drange = dt * c/(sqrt(er_ice)*2);
      elev_uniform = (elev_max:-drange:elev_min).';
      % update image_data
      Nt = size(img,1);
      img = [img;...
        zeros(length(elev_uniform)-Nt,size(img,2))];
      warning('off','MATLAB:interp1:NaNinY')
      for idx = 1:length(surface)
        range = min(time,surface(idx)) * c/2 ...
          + max(0,time-surface(idx)) * c/(sqrt(er_ice)*2);
        elev = elevation(idx) - range;
        img(:,idx) = interp1(elev,...
          img(1:Nt,idx),elev_uniform,'linear');
      end
      warning('on','MATLAB:interp1:NaNinY')
      pass(pass_idx).data_elev = img;
      pass(pass_idx).data_elev_yaxis = elev_uniform;
      imagesc(pass(master_idx).along_track/1e3, elev_uniform, img);
      ylabel('Elevation (m, WGS-84)');
      xlabel('Along-track (km)');
      set(gca,'YDir','normal');
      
      % Convert layerPnts_y from TWTT to WGS_84 Elevation
      elevation = pass(master_idx).elev;
      surface = pass(pass_idx).layers(1).twtt_ref;
      for lay_idx = 1:length(pass(pass_idx).layers)
        layer_y_curUnit = pass(pass_idx).layers(lay_idx).twtt_ref;
        for pnt_idx = 1:length(layer_y_curUnit)
          range = min(layer_y_curUnit(pnt_idx),surface(pnt_idx))*c/2 ...
            +max(0,layer_y_curUnit(pnt_idx)-surface(pnt_idx)) * c/(sqrt(er_ice)*2);
          layer_y_curUnit(pnt_idx) = elevation(pnt_idx) - range;
        end
        layer_y_curUnit(isnan(pass(pass_idx).layers(lay_idx).twtt_ref)) = NaN;
        pass(pass_idx).layers(lay_idx).layer_elev = layer_y_curUnit;
        hold on;
        plot(pass(master_idx).along_track/1e3, pass(pass_idx).layers(lay_idx).layer_elev);
      end
  
    else
      imagesc(lp(data(rbins,:,pass_idx)))
      ylabel('Range bin');
      xlabel('Range line');
    end
    colormap(1-gray(256));
    title(sprintf('%s_%03d %d',pass(pass_idx).param_sar.day_seg,pass(pass_idx).param_multipass.cmd.frms(1),pass(pass_idx).direction),'interpreter','none')
    %caxis([-90 8]);
  else
    % Form interferogram (couple options)
    complex_data = fir_dec(data(rbins,:,pass_idx) .* conj(data(rbins,:,master_idx)),ones(1,11)/11,1);
    if ~exist('equalization_rlines','var') || isempty(equalization_rlines)
      new_equalization(pass_idx) = mean(complex_data(:)); % equalization only valid when motion compensation with phase is used
    else
      new_equalization(pass_idx) = mean(mean(complex_data(:,equalization_rlines))); % equalization only valid when motion compensation with phase is used
    end
    if param.multipass.comp_mode == 1
      complex_data = complex_data ./ new_equalization(pass_idx);
    end
    % Plot interferogram
    if param.multipass.comp_mode == 4
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
linkaxes(h_data_axes,'xy');

%% Save Result
% =========================================================================
[fn_dir,fn_name] = fileparts(fn);
fn_multipass = fullfile(fn_dir, sprintf('%s_multipass%02.0f.mat',fn_name, param.multipass.comp_mode));
param_sar = pass(master_idx).param_sar;
param_records = pass(master_idx).param_records;

if ~isfield(pass,'param_multipass')
  param_multipass = read_param_xls(param_fn,param_combine_passes.combine_passes.passes(master_idx).day_seg); %reads parameter sheet for given pass
else
  param_multipass = pass.param_multipass; %reads parameter sheet for given pass
end
param_multipass = merge_structs(param,param_multipass);

fprintf('Saving %s (%s)\n', fn_multipass, datestr(now));
save(fn_multipass,'-v7.3','pass','data','ref','param_records','param_sar','param_multipass');

if param.multipass.comp_mode ~= 2 || strcmp('echo',param.multipass.echo_sar)
  return
end

%% Array Processing (comp_mode == 2)
% Package data to call array_proc.m
% 1. Data
% 2. Trajectory and attitude
% 3. Array processing parameters
data = {permute(data,[1 2 4 5 3])};

param.array = [];
param.array.method = 1;
% param.array.Nsv = 256; param.array = rmfield(param.array,'theta'); % Forces default theta
% param.array.theta = linspace(-20,20,256);
param.array.theta = linspace(-6,6,256);
param.array.Nsrc = 2;
param.array.bin_rng = [-4:4];
param.array.line_rng = [-20:20];
param.array.dbin = 1;
param.array.dline = 11;
param.array.freq_rng = 1;
h_fig_baseline = figure(200); clf;
h_plot_baseline = [];
h_legend_baseline = {};
wf_cells = {[],[],[]};
for pass_out_idx = 1:length(pass_en_idxs)
  pass_idx = pass_en_idxs(pass_out_idx);
  
  for wf_id = 1:length(pass(pass_out_idx).wfs)
    wf_cells{wf_id} = [wf_cells{wf_id}; ...
      wf_id, pass_out_idx];
  end
  
  param.array.fcs{1}{pass_out_idx}.pos = along_track;
  param.array.fcs{1}{pass_out_idx}.pos(2,:) = pass(pass_idx).ref_y;
  param.array.fcs{1}{pass_out_idx}.pos(3,:) = pass(pass_idx).ref_z;
  param.array.fcs{1}{pass_out_idx}.base_line ...        
    = sqrt( (pass(pass_idx).ref_z - pass(baseline_master_idx).ref_z).^2 ...
      + (pass(pass_idx).ref_y - pass(baseline_master_idx).ref_y).^2 );
    
  h_plot_baseline(end+1) = plot(param.array.fcs{1}{pass_out_idx}.base_line);
  h_legend_baseline{end+1} = sprintf('%d',pass_idx);
  hold on;

  param.array.fcs{1}{pass_out_idx}.surface = ref.surface;
end
xlabel('Range line');
ylabel('Baseline (m)');
grid on;
legend(h_plot_baseline,h_legend_baseline);

param.array.wfs.time = ref.time;
dt = param.array.wfs.time(2)-param.array.wfs.time(1);
param.array_proc.bin0 = param.array.wfs.time(1)/dt;
param.array.sv_fh = @array_proc_sv;
param.array.wfs.fc = ref.wfs(ref.wf).fc;
param.array.imgs = {[ones(length(pass_en_idxs),1), (1:length(pass_en_idxs)).']};
% param.array.imgs = wf_cells;
param.array.tomo_en = true;

%%
array_proc_methods;
param = array_proc(param);
param.array.method = STANDARD_METHOD;
[param_array0,result0] = array_proc(param,{data{1}(:,:,:,:,pass_en_idxs)});
% param.array.method = MVDR_METHOD;
% [param_array1,result1] = array_proc(param,{data{1}(:,:,:,:,pass_en_idxs)});
param.array.method = MUSIC_METHOD;
[param_array2,result2] = array_proc(param,{data{1}(:,:,:,:,pass_en_idxs)});


%% Save Results
[fn_dir,fn_name] = fileparts(fn);

Tomo = result0.tomo;
Data = result0.img;
GPS_time = ref.gps_time(param_array0.array_proc.lines);
Latitude = ref.lat(param_array0.array_proc.lines);
Longitude = ref.lon(param_array0.array_proc.lines);
Elevation = ref.elev(param_array0.array_proc.lines);
Roll = ref.roll(param_array0.array_proc.lines);
Pitch = ref.pitch(param_array0.array_proc.lines);
Heading = ref.heading(param_array0.array_proc.lines);
Surface = ref.surface(param_array0.array_proc.lines);
Bottom = nan(size(Surface));
param_sar = pass(baseline_master_idx).param_sar;
param_records = pass(baseline_master_idx).param_records;
param_array = param_array0;
Time = pass(master_idx).time(param_array2.array_proc.bins);
file_version = '1';
fn_mat = fullfile(fn_dir,[fn_name output_fn_midfix '_standard.mat']);
save('-v7.3',fn_mat,'Tomo','Data','Latitude','Longitude','Elevation','GPS_time', ...
  'Surface','Bottom','Time','param_array','param_records', ...
  'param_sar', 'Roll', 'Pitch', 'Heading', 'file_version');

Tomo = result2.tomo;
Data = result2.img;
GPS_time = ref.gps_time(param_array2.array_proc.lines);
Latitude = ref.lat(param_array2.array_proc.lines);
Longitude = ref.lon(param_array2.array_proc.lines);
Elevation = ref.elev(param_array2.array_proc.lines);
Roll = ref.roll(param_array2.array_proc.lines);
Pitch = ref.pitch(param_array2.array_proc.lines);
Heading = ref.heading(param_array2.array_proc.lines);
Surface = ref.surface(param_array2.array_proc.lines);
Bottom = nan(size(Surface));
param_sar = pass(baseline_master_idx).param_sar;
param_records = pass(baseline_master_idx).param_records;
param_array = param_array2;
Time = pass(master_idx).time(param_array2.array_proc.bins);
file_version = '1';
fn_mat = fullfile(fn_dir,[fn_name output_fn_midfix '_music.mat']);
save('-v7.3',fn_mat,'Tomo','Data','Latitude','Longitude','Elevation','GPS_time', ...
  'Surface','Bottom','Time','param_array','param_records', ...
  'param_sar', 'Roll', 'Pitch', 'Heading', 'file_version');



return

