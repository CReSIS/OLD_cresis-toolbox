function echo_info = publish_echogram(param,mdata,lay)
% echo_info = publish_echogram(param,mdata,lay)
%
% param: structure controlling the plotting of the echogram
%  .fig_hand: figure handle to use for plot
%  .num_x_tics: number of xtics across the bottom (the xtics give
%     the geographic info), usually 4 for figure copy and 6 for
%     full page print out
%  .depth: 1x2 double vector representing the limits of the depth
%     axis to plots (-500 to 3500 typical, deep antarctica -500 to 4000,
%     but for custom plots -100 to 2000 or whatever highlights the data
%     best
%  .elev_comp: Elevation compensation mode
%     0: no elevation compensation (twtt)
%     1: level flight elevation compensation
%     2: depth elevation compensation
%     3: WGS-84 elevation compensation
%  .er_ice: the permittivity used to create the depth axis, default is 
%     loaded from physical_constants.
%  .plot_quality: boolean scalar, default is false, if true, plots each
%     layer with colors corresponding to the layers quality
% echo_info: structure with handles to plots/axes
%  .ah_echo_time: time axis handle
%  .ah_echo: depth axis handle
%  .h_surf: surface plot handle
%  .h_bot: bottom plot handle
%
% Example: see run_publish_echogram
%
% Author: John Paden
%
% See also: publish_map, run_publish_echogram

physical_constants;
if ~isfield(param,'er_ice') || isempty(param.er_ice)
  param.er_ice = er_ice;
end
if ~isfield(param,'er_depth') || isempty(param.er_depth)
  param.er_depth = 0;
end
if ~isfield(param,'axis_type') || isempty(param.axis_type)
  param.axis_type = 'standard';
end
if ~isfield(param,'plot_quality') || isempty(param.plot_quality) 
  param.plot_quality = false;
end
if ~isfield(param,'caxis') || isempty(param.caxis) 
  param.caxis = [];
end
if ~isfield(param,'depth_offset') || isempty(param.depth_offset)
  % Adds an offset to the depth/elevation/range tics (in meters)
  param.depth_offset = 0;
end
if ~isfield(param,'time_offset') || isempty(param.time_offset) 
  % Adds an offset to the two way travel time tics (in seconds)
  param.time_offset = 0;
end
if ~isreal(mdata.Data)
  warning('Input data are complex. Taking the abs()^2 of the data.');
  mdata.Data = abs(mdata.Data).^2;
end

% =======================================================================
% Create Surface and Bottom Variables
% =======================================================================
lay.Bottom  = lay.layerData{2}.value{2}.data;
lay.Surface = lay.layerData{1}.value{2}.data;
lay.Thickness = lay.Bottom-lay.Surface;
if param.plot_quality
  lay.Surface_Quality = lay.layerData{1}.quality;
  lay.Bottom_Quality = lay.layerData{2}.quality;
end
neg_idxs = find(lay.Thickness < 0 & isfinite(lay.Thickness));
if ~isempty(neg_idxs)
  warning('  Negative thickness detected');
  lay.Bottom(neg_idxs) = lay.Surface(neg_idxs);
end

% Interpolate layer onto GPS times of echogram
warning('off','MATLAB:interp1:NaNinY');
if param.plot_quality
  lay.Surface_Quality = interp1(lay.GPS_time,lay.Surface_Quality,mdata.GPS_time,'linear','extrap');
  lay.Bottom_Quality = interp1(lay.GPS_time,lay.Bottom_Quality,mdata.GPS_time,'linear','extrap');
end
lay.Surface = interp1(lay.GPS_time,lay.Surface,mdata.GPS_time,'linear','extrap');
lay.Bottom = interp1(lay.GPS_time,lay.Bottom,mdata.GPS_time,'linear','extrap');
warning('on','MATLAB:interp1:NaNinY');

% Create fast-time correction vector
%   The layer file may contain a different elevation profile
%   than the data file. The surface and bottom need to be adjusted
%   to account for these differences.
elev_interp    = interp1(lay.GPS_time,lay.Elevation,mdata.GPS_time,'linear','extrap');
fast_time_correction = (mdata.Elevation - elev_interp)/(c/2);

% =======================================================================
% Elevation compensation
% =======================================================================

if param.elev_comp == 1
  %% Relative elevation compensation
  max_elev = max(mdata.Elevation);
  dRange = max_elev - mdata.Elevation;
  dt = mdata.Time(2)-mdata.Time(1);
  dBins = round(dRange / (c/2) / dt);
  dtime = dRange/(c/2);
  zero_pad_len = max(abs(dBins));
  mdata.Data = cat(1,mdata.Data,zeros(zero_pad_len,size(mdata.Data,2)));
  mdata.Time = mdata.Time(1) + (mdata.Time(2)-mdata.Time(1)) * (0:1:size(mdata.Data,1)-1);

  warning off
  for rline = 1:size(mdata.Data,2)   
    mdata.Data(:,rline) = interp1(mdata.Time, mdata.Data(:,rline), mdata.Time - dtime(rline), 'linear',0);
    mdata.Elevation(rline) = mdata.Elevation(rline) + dRange(rline);
    lay.Surface(rline) = lay.Surface(rline) + dtime(rline);
    lay.Bottom(rline) = lay.Bottom(rline) + dtime(rline);
  end
  warning on;
elseif param.elev_comp == 2
  %% Compensate for surface variations (i.e. flatten to low pass filtered
  % version of surface)
  
  % Filter out high frequencies from surface layer
  lay.Surface_Filled = interp_finite(lay.Surface,0); % Need surface points everywhere for filtering operation
  if ~isfield(param,'surf_filt_en') || isempty(param.surf_filt_en) || param.surf_filt_en
    if length(lay.Surface_Filled) >= 100
      [Bfilt,Afilt] = butter(2,0.02);
      surf_filt = filtfilt(Bfilt,Afilt,lay.Surface_Filled);
      tmp = polyval(polyfit(1:51,reshape(lay.Surface_Filled(1:51),[1 51]),2),1:51);
      surf_filt(1:50) = tmp(1:50) - tmp(51) + surf_filt(51);
      tmp = polyval(polyfit(1:51,reshape(lay.Surface_Filled(end-50:end),[1 51]),2),1:51);
      surf_filt(end-49:end) = tmp(2:51) - tmp(1) + surf_filt(end-50);
    else
      surf_filt = medfilt1(lay.Surface_Filled,round(length(lay.Surface_Filled)/20)*2+1);
    end
  elseif ~param.surf_filt_en || isempty(param.surf_filt_en)
    surf_filt = lay.Surface_Filled;
  end

  lay.Surface = surf_filt;
  mdata.Elevation = mean(mdata.Elevation) + (surf_filt - mean(surf_filt)) * c/2;  
  

  Surface_Depth = zeros(size(mdata.GPS_time));
  depth_range = eval(param.depth);

  
    % Determine sampling
  max_er = max(param.er_ice);
  dt = mdata.Time(2) - mdata.Time(1);
  dz = dt * c/2/sqrt(max_er);

  % Depth axis
  depth = (round(depth_range(1)/dz)*dz : dz : depth_range(end)).';
  
  % Depth time axis
  depth_time = size(depth);
  % Above surface
  depth_time(depth <= 0) = depth(depth <= 0) / (c/2);
  % Below surface and within dielectric profile defined depth
  if length(param.er_depth) > 1
    TWtime = genPropProfileFromPerm(param.er_depth,param.er_ice,1);
    profile_idxs = depth > 0 & depth < param.er_depth(end);
    depth_time(profile_idxs) = interp1(param.er_depth, [0; TWtime], depth(profile_idxs));
  else
    TWtime = 0;
  end
  % Below surface and below dielectric profile defined depth
  const_idxs = depth >= param.er_depth(end);
  depth_time(const_idxs) = TWtime(end) + (depth(const_idxs) - param.er_depth(end)) / (c/2/sqrt(param.er_ice(end)));
  if length(depth_time) > length(mdata.Time);
      mdata.Data = cat(1,mdata.Data,ones(length(depth_time) -length(mdata.Time),size(mdata.Data,2))*mean(mdata.Data(end,:)));
      mdata.Time = cat(1,mdata.Time,mdata.Time(end) + (depth_time(end)-depth_time(end-1))*(1:(length(depth_time) -length(mdata.Time)))');
  end
  % Re-interpolate data to a constant depth axis
  newData = zeros(length(depth),size(mdata.Data,2));
  warning off;
  mdata.Data(mdata.Data == 0) = NaN;
  for rline = 1:size(mdata.Data,2)
    newData(:,rline) = interp1(mdata.Time, mdata.Data(:,rline), ...
      lay.Surface(rline) + depth_time);
    lay.Bottom(:) = lay.Bottom(rline) - lay.Surface(rline);
  end
  mdata.Data = newData;
  mdata.Data(isnan(mdata.Data)) = 0;
  warning on;

elseif param.elev_comp == 3
  %% Elevation compensate to WGS-84 y-axis
  
  % Filter surface
  lay.Surface_Filled = interp_finite(lay.Surface,0); % Need surface points everywhere for filtering operation
  if ~isfield(param,'surf_filt_en') || param.surf_filt_en
    if length(lay.Surface_Filled) >= 100
      [Bfilt,Afilt] = butter(2,0.02);
      surf_filt = filtfilt(Bfilt,Afilt,lay.Surface_Filled);
      tmp = polyval(polyfit(1:51,reshape(lay.Surface_Filled(1:51),[1 51]),2),1:51);
      surf_filt(1:50) = tmp(1:50) - tmp(51) + surf_filt(51);
      tmp = polyval(polyfit(1:51,reshape(lay.Surface_Filled(end-50:end),[1 51]),2),1:51);
      surf_filt(end-49:end) = tmp(2:51) - tmp(1) + surf_filt(end-50);
    else
      surf_filt = medfilt1(lay.Surface_Filled,round(length(lay.Surface_Filled)/20)*2+1);
    end
  elseif ~param.surf_filt_en || isempty(param.surf_filt_en)
    surf_filt = lay.Surface_Filled;
  end
  
  % Remove data before zero time
  negative_bins = mdata.Time < 0;
  mdata.Time = mdata.Time(~negative_bins);
  mdata.Data = mdata.Data(~negative_bins,:);
  % Create elevation axis to interpolate to
  max_elev = max(mdata.Elevation);
  min_elev = min(mdata.Elevation - surf_filt*c/2 - (mdata.Time(end)-surf_filt)*c/2/sqrt(param.er_ice));
  dt = mdata.Time(2)-mdata.Time(1);
  dr = dt * c/2 / sqrt(param.er_ice);
  elev_axis = max_elev:-dr:min_elev;

  % Zero pad data to create space for interpolated data
  zero_pad_len = length(elev_axis) - length(mdata.Time);
  mdata.Data = cat(1,mdata.Data,zeros(zero_pad_len,size(mdata.Data,2)));

  % Determine the corrections to apply to elevation and layers
  dRange = max_elev - mdata.Elevation;
  dBins = round(dRange / (c/2) / dt);
  dtime = dRange/(c/2);

  warning off
  for rline = 1:size(mdata.Data,2)
    % Determine elevation bins before surface
    surf_elev = mdata.Elevation(rline) - surf_filt(rline) * c/2;
    dt_air = dr/(c/2);
    time0 = -(max_elev - mdata.Elevation(rline))/(c/2);
    last_air_idx = find(elev_axis > surf_elev,1,'last');
    new_time = (time0 + dt_air*(0:last_air_idx-1)).';
    if isempty(last_air_idx)
      % Radar is on the surface of the ice
      last_air_idx = 0;
    end
    if last_air_idx < length(elev_axis)
      % Determine elevation bins after surface
      dt_ice = dr/(c/2/sqrt(param.er_ice));
      first_ice_idx = last_air_idx + 1;
      time0 = surf_filt(rline) + (surf_elev - elev_axis(first_ice_idx))/(c/2/sqrt(param.er_ice));
      new_time = cat(1,new_time, (time0 + dt_ice*(0:length(elev_axis)-length(new_time)-1)).');
    end
    mdata.Data(:,rline) = interp1(mdata.Time, mdata.Data(1:length(mdata.Time),rline), new_time, 'linear',0);
    mdata.Elevation(rline) = mdata.Elevation(rline) + dRange(rline);
    lay.Surface(rline) = lay.Surface(rline) + dtime(rline);
    lay.Bottom(rline) = lay.Bottom(rline) + dtime(rline);
  end
  warning on

end

% Create depth axis
good_surface_vals = lay.Surface + fast_time_correction;
good_surface_vals = good_surface_vals(isfinite(good_surface_vals));
mean_surface_time = mean(good_surface_vals);
if ~isfinite(mean_surface_time)
  mean_surface_time = 0;
end

% Limit depths according to input param.depth
if param.elev_comp == 3
  DSurface = mdata.Elevation - lay.Surface*c/2;
  DBottom = mdata.Elevation - lay.Surface*c/2 - (lay.Bottom-lay.Surface)*c/2/sqrt(param.er_ice);
  Surface_Elev = DSurface;
  Bbad = sum(isnan(DBottom)) / numel(DBottom);
  % Example: param.depth = '[min(Surface_Elev) - 15 max(Surface_Elev)+3]';
  % Example: param.depth = '[100 120]';
  % Example: param.depth = '[publish_echogram_switch(Bbad,0.25,Surface_Elev,-1600,DBottom,-100),max(Surface_Elev+50)]';

  depth_range = eval(param.depth);
  depth_good_idxs = find(elev_axis >= depth_range(1) & elev_axis <= depth_range(end));
  if isfield(param,'detrend') && ~isempty(param.detrend) && strcmpi(param.detrend.mode,'polynomial')
    detrend_depth_range = eval(param.detrend.depth);
    detrend_depth_good_idxs = find(elev_axis >= detrend_depth_range(1) ...
      & elev_axis <= detrend_depth_range(end));
  end
elseif param.elev_comp == 2
  Depth = depth;
  DSurface = lay.Surface; % All zero
  DBottom = interp1(depth_time,depth,lay.Bottom);
  Surface_Depth = DSurface;
  Bbad = sum(isnan(DBottom)) / numel(DBottom);
  
  depth_good_idxs = 1:length(Depth);
  if isfield(param,'detrend') && ~isempty(param.detrend) && strcmpi(param.detrend.mode,'polynomial')
    detrend_depth_range = eval(param.detrend.depth);
    detrend_depth_good_idxs = find(Depth >= detrend_depth_range(1) ...
      & Depth <= detrend_depth_range(end));
  end
else
  Depth = (mdata.Time-mean_surface_time)*c/2/sqrt(param.er_ice);
  
  DBottom = lay.Bottom - mean_surface_time + fast_time_correction;
  DBottom = DBottom*c/2/sqrt(param.er_ice);
  DSurface = lay.Surface - mean_surface_time + fast_time_correction;
  DSurface = DSurface*c/2/sqrt(param.er_ice);
  Surface_Depth = DSurface;
  Bbad = sum(isnan(DBottom)) / numel(DBottom);
  
  depth_range = eval(param.depth);
  depth_good_idxs = find(Depth >= depth_range(1) & Depth <= depth_range(end));
  if isfield(param,'detrend') && ~isempty(param.detrend) && strcmpi(param.detrend.mode,'polynomial')
    detrend_depth_range = eval(param.detrend.depth);
    detrend_depth_good_idxs = find(Depth >= detrend_depth_range(1) ...
      & Depth <= detrend_depth_range(end));
  end
end

% =======================================================================
% Image filtering
% =======================================================================
if isfield(param,'filter') && ~isempty(param.filter)
  %% Apply optional image filter
  mdata.Data = param.filter(mdata.Data);
end

% =======================================================================
% Detrend
% =======================================================================
if isfield(param,'detrend') && ~isempty(param.detrend)
  if strcmpi(param.detrend.mode,'tonemap')
    %% Automated detrending Code using tonemap
    
    hdr1(:,:,1) = abs(mdata.Data(depth_good_idxs,:));
    hdr1(:,:,2) = abs(mdata.Data(depth_good_idxs,:));
    hdr1(:,:,3) = abs(mdata.Data(depth_good_idxs,:));
    % for generating synthetic HDR images
    detrend_tonemap = tonemap(hdr1, 'AdjustLightness', [0.1 1], 'AdjustSaturation', 1.5);
    detrend_tonemap = detrend_tonemap(:,:,2);

  elseif strcmpi(param.detrend.mode,'polynomial')
    %% Polynomial Detrending Code
    % Detrending splits echogram into 3 sections. The assumption is that
    % elevation compensation is set to 2 (i.e. surface is flattened)
    % Three sections are made to be continuous:
    % 1. Constant (usually above the ice surface)
    % 2. Polynomial fit to log of mean range bin power
    % 3. Constant (usually the noise floor)
    
    if param.elev_comp ~= 2
      warning('Currently polynomial detrending assumes elevation compensation set to 2 (surface flattening)');
    end
    
    detrend_rbins = detrend_depth_good_idxs; % Specifies bins for the polynomial section
    detrend_guard = round(0.1*length(detrend_depth_good_idxs));
    detrend_poly_order = param.detrend.poly_order; % Specify order of second section
    
    detrend_rlines = 1:size(mdata.Data,2);
    
    % Get the log of the mean range bin power
    detrend_profile = lp(mean(abs(mdata.Data),2));
    
    % Create an xaxis to prevent badly conditioned polynomials
    if isempty(detrend_rbins)
      detrend_rbins = (1:length(detrend_profile));
    end
    xaxis = (detrend_rbins - mean(detrend_rbins));
    
    % Find the polynomial coefficients
    detrend_poly = polyfit(xaxis, reshape(detrend_profile(detrend_rbins),size(xaxis)), detrend_poly_order);
    
    % Section 2: Evaluate the polynomial
    detrend_profile(detrend_rbins) = polyval(detrend_poly,xaxis);
    
    % Section 1: Constant
    detrend_profile(1:detrend_rbins(1)-1) = detrend_profile(detrend_rbins(1));
    
    % Section 3: Constant
    detrend_profile(detrend_rbins(end)-detrend_guard:end) = detrend_profile(detrend_rbins(end)-detrend_guard);
    
    % Check detrending
    if 0
      figure;
      plot(lp(mean(abs(mdata.Data),2)),'r');
      hold on;
      plot(detrend_profile);
      hold off;
    end
    
    % Apply detrending
    if all(isfinite(detrend_profile))
      mdata.Data = mdata.Data ./ repmat(10.^(detrend_profile/10), [1 size(mdata.Data,2)]);
    end
   end
end

% Create echogram plot
figure(param.fig_hand); clf;

ah_echo_time = axes; % 2-way travel time axis
ah_echo = axes; % Depth axis with data
if param.elev_comp == 3
  %% WGS-84 Elevation elevation comp plotting
  echogram_vals = lp(mdata.Data(depth_good_idxs,:));
  
  if isfield(param,'detrend') && ~isempty(param.detrend) && strcmpi(param.detrend.mode,'tonemap')
    echo_info.image = imagesc([],elev_axis(depth_good_idxs)+param.depth_offset, ...
      detrend_tonemap,'Parent',ah_echo);
    axis(ah_echo_time,[0.5 size(detrend_tonemap,2)+0.5 reshape(new_time(depth_good_idxs([1 end]))*1e6 + param.time_offset*1e6,[1 2])])
  else
    echo_info.image = imagesc([],elev_axis(depth_good_idxs)+param.depth_offset, ...
      echogram_vals,'Parent',ah_echo);
    axis(ah_echo_time,[0.5 size(echogram_vals,2)+0.5 reshape(new_time(depth_good_idxs([1 end]))*1e6 + param.time_offset*1e6,[1 2])])
  end
  set(ah_echo,'YDir','Normal');
  ylabel(ah_echo,sprintf('WGS-84 Elevation, e_r = %.2f (m)', param.er_ice));
  
else
  %% Everything except WGS-84 Elevation elevation comp plotting
  echogram_vals = lp(mdata.Data(depth_good_idxs,:));
  
  if isfield(param,'detrend') && ~isempty(param.detrend) && strcmpi(param.detrend.mode,'tonemap')
    echo_info.image = imagesc([],Depth(depth_good_idxs)+param.depth_offset, ...
      detrend_tonemap,'Parent',ah_echo);
    axis(ah_echo_time,[0.5 size(detrend_tonemap,2)+0.5 reshape(mdata.Time(depth_good_idxs([1 end]))*1e6 + param.time_offset*1e6,[1 2])])
  else
    echo_info.image = imagesc([],Depth(depth_good_idxs)+param.depth_offset, ...
      echogram_vals,'Parent',ah_echo);
    axis(ah_echo_time,[0.5 size(echogram_vals,2)+0.5 reshape(mdata.Time(depth_good_idxs([1 end]))*1e6 + param.time_offset*1e6,[1 2])])
  end
  if length(param.er_ice) == 1;
    ylabel(ah_echo,sprintf('depth, e_r = %.2f (m)', param.er_ice));
  else
    ylabel(ah_echo,sprintf('depth, e_r from profile (m)', param.er_ice));
  end
end
ylabel(ah_echo_time,'Propagation delay (us)');
set(ah_echo_time,'XTickLabel','');
set(ah_echo_time,'YAxisLocation','Right');
set(ah_echo_time,'YDir','Reverse');
if strcmpi(param.axis_type,'bars')
  xlabel(ah_echo,'');
  ylabel(ah_echo,'');
  set(ah_echo,'XTick',[])
  set(ah_echo,'YTick',[])
  set(ah_echo_time,'Visible','off');
  xlims = xlim;
  ylims = ylim;
  
  ybar_length_exact = diff(ylims)*0.2;
  [~,min_log_match] = min(abs(mod(log(ybar_length_exact),log(10)) - log(1:10)));
  ybar_length = min_log_match * 10^floor(log(ybar_length_exact) / log(10));
  ybar_length_norm = ybar_length/diff(ylims);

  along_track = geodetic_to_along_track(mdata.Latitude,mdata.Longitude,mdata.Elevation);
  xbar_length_exact = (along_track(end) - along_track(1))*0.2;
  [~,min_log_match] = min(abs(mod(log(xbar_length_exact),log(10)) - log(1:10)));
  xbar_length = min_log_match * 10^floor(log(xbar_length_exact) / log(10));
  xbar_length_norm = xbar_length/(along_track(end) - along_track(1));
  
  set(ah_echo,'Units','pixels')
  pixel_pos = get(ah_echo,'Position');
  set(ah_echo,'Units','normalized')
  xbar_width = 4 / pixel_pos(4);
  ybar_width = 4 / pixel_pos(3);
  
  if strcmpi(get(ah_echo,'YDir'),'reverse')
    echo_info.yaxis_bar = patch(xlims(1)+diff(xlims)*ybar_width*[2 3 3 2 2], ylims(2)-diff(ylims)*(xbar_width*4+[0 0 ybar_length_norm ybar_length_norm 0]),'k','Parent',ah_echo);
    echo_info.xaxis_bar = patch(xlims(1)+diff(xlims)*(ybar_width*4+[0 xbar_length_norm xbar_length_norm 0 0]),ylims(2)-diff(ylims)*xbar_width*[2 2 3 3 2],'k','Parent',ah_echo);
    echo_info.yaxis_bar_text = text(xlims(1)+diff(xlims)*ybar_width*4, ...
      ylims(2)-diff(ylims)*(xbar_width*4+ybar_length_norm/2), ...
      sprintf('%.0f m',ybar_length),'Parent',ah_echo,'Color',[0 0 0],'FontWeight','Bold');
    if xbar_length >= 1000
      echo_info.xaxis_bar_text = text(xlims(1)+diff(xlims)*(ybar_width*4+xbar_length_norm/2), ...
        ylims(2)-diff(ylims)*xbar_width*3, ...
        sprintf('%.0f km',xbar_length/1000),'Parent',ah_echo,'Color',[0 0 0], ...
        'HorizontalAlignment','Center','VerticalAlignment','Bottom','FontWeight','Bold');
    else
      echo_info.xaxis_bar_text = text(xlims(1)+diff(xlims)*(ybar_width*4+xbar_length_norm/2), ...
        ylims(2)-diff(ylims)*xbar_width*3, ...
        sprintf('%.0f m',xbar_length),'Parent',ah_echo,'Color',[0 0 0], ...
        'HorizontalAlignment','Center','VerticalAlignment','Bottom','FontWeight','Bold');
    end
  else
    echo_info.yaxis_bar = patch(xlims(1)+diff(xlims)*ybar_width*[2 3 3 2 2], ylims(1)+diff(ylims)*(xbar_width*4+[0 0 ybar_length_norm ybar_length_norm 0]),'k','Parent',ah_echo);
    echo_info.xaxis_bar = patch(xlims(1)+diff(xlims)*(ybar_width*4+[0 xbar_length_norm xbar_length_norm 0 0]),ylims(1)+diff(ylims)*xbar_width*[2 2 3 3 2],'k','Parent',ah_echo);
    echo_info.yaxis_bar_text = text(xlims(1)+diff(xlims)*ybar_width*4, ...
      ylims(1)+diff(ylims)*(xbar_width*4+ybar_length_norm/2), ...
      sprintf('%.0f m',ybar_length),'Parent',ah_echo,'Color',[0 0 0],'FontWeight','Bold');
    if xbar_length >= 1000
      echo_info.xaxis_bar_text = text(xlims(1)+diff(xlims)*(ybar_width*4+xbar_length_norm/2), ...
        ylims(1)+diff(ylims)*xbar_width*3, ...
        sprintf('%.0f km',xbar_length/1000),'Parent',ah_echo,'Color',[0 0 0], ...
        'HorizontalAlignment','Center','VerticalAlignment','Bottom','FontWeight','Bold');
    else
      echo_info.xaxis_bar_text = text(xlims(1)+diff(xlims)*(ybar_width*4+xbar_length_norm/2), ...
        ylims(1)+diff(ylims)*xbar_width*3, ...
        sprintf('%.0f m',xbar_length),'Parent',ah_echo,'Color',[0 0 0], ...
        'HorizontalAlignment','Center','VerticalAlignment','Bottom','FontWeight','Bold');
    end
  end
else
  xtl = create_standard_x_labels(mdata.Latitude,mdata.Longitude,mdata.Elevation,param.num_x_tics);
  add_x_labels(ah_echo,xtl,{'dist','lat','lon'});
  set(ah_echo_time,'Units','normalized');
  set(ah_echo,'Units','normalized');
  set(ah_echo_time,'Position',get(ah_echo,'Position'));
end

if isfield(param,'colormap') && ~isempty(param.colormap)
  colormap(param.colormap);
else
  colormap(1-gray);
end
if ~isempty(param.caxis)
  if isfield(param,'detrend') && ~isempty(param.detrend) && strcmpi(param.detrend.mode,'polynomial')
    check_vals = lp(mdata.Data(detrend_depth_good_idxs,:));
    check_vals = check_vals(isfinite(check_vals));
    check_vals = sort(check_vals(:));
    if length(check_vals) >= 2
      new_caxis = check_vals(1+round(param.caxis*(length(check_vals)-1)));
      caxis(new_caxis);
    end
  else
    echogram_vals = sort(echogram_vals(:));
    echogram_vals = echogram_vals(isfinite(echogram_vals));
    if length(echogram_vals) >= 2
      new_caxis = echogram_vals(1+round(param.caxis*(length(echogram_vals)-1)));
      caxis(new_caxis);
    end
  end
end

echo_info.h_title = title(ah_echo,sprintf('Data Frame ID: %s', param.frm_id),'Interpreter','none','FontWeight','normal');

% Plot surface and bottom layer
hold on
if param.plot_quality
  
  tmp_layer = DSurface;
  tmp_layer(lay.Surface_Quality~=1) = NaN;
  echo_info.h_surf(1) = plot(tmp_layer,'g--');
  
  tmp_layer = DSurface;
  tmp_layer(lay.Surface_Quality~=2) = NaN;
  echo_info.h_surf(2) = plot(tmp_layer,'y--');
  
  tmp_layer = DSurface;
  tmp_layer(lay.Surface_Quality~=3) = NaN;
  echo_info.h_surf(3) = plot(tmp_layer,'r--');
  
  tmp_layer = DBottom;
  tmp_layer(lay.Bottom_Quality~=1) = NaN;
  echo_info.h_bot(1) = plot(tmp_layer,'g--');
  
  tmp_layer = DBottom;
  tmp_layer(lay.Bottom_Quality~=2) = NaN;
  echo_info.h_bot(2) = plot(tmp_layer,'y--');
  
  tmp_layer = DBottom;
  tmp_layer(lay.Bottom_Quality~=3) = NaN;
  echo_info.h_bot(3) = plot(tmp_layer,'r--');
else
  echo_info.h_surf = plot(ah_echo,DSurface,'--m');
  echo_info.h_bot = plot(ah_echo,DBottom,'--r');
end
hold off

echo_info.ah_echo_time = ah_echo_time;
echo_info.ah_echo = ah_echo;

end

function depth_range = publish_echogram_switch(Bbad,Bbad_Threshold,Surface_Elev,Surface_Elev_Offset,DBottom,DBottom_Offset)

if Bbad > Bbad_Threshold
  depth_range = min(Surface_Elev+Surface_Elev_Offset);
else
  depth_range = min(DBottom+DBottom_Offset);
end

end
