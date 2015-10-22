if 0
  
  tic;
  physical_constants;
  
  % =====================================================================
  % User Settings
  % =====================================================================
  
  fns = {};
  new_freq_rng = {};
  swath_mode = {};
  
  base_path = '/mnt/products/rds/2013_Antarctica_Basler/CSARP_music2/20131220_03';
  base_path2 = '/mnt/products/rds/2013_Antarctica_Basler/CSARP_tomo/20131220_03';
  for frm = [3 4 5 6]
    fns{end+1} = fullfile(base_path,sprintf('Data_img_01_20131220_03_%03i',frm));
    new_freq_rng{end+1} = [1:64];
    swath_mode{end+1} = 2;
    fns{end+1} = fullfile(base_path,sprintf('Data_img_02_20131220_03_%03i',frm));
    new_freq_rng{end+1} = [33:64];
    swath_mode{end+1} = 1;
    fns{end+1} = fullfile(base_path,sprintf('Data_img_03_20131220_03_%03i',frm));
    new_freq_rng{end+1} = [1:33];
    swath_mode{end+1} = 1;
  end
  
  for fn_idx = 1:length(fns)
    fn = fns{fn_idx};
    
    fprintf('Reading %s\n', fn);
    %% Load 3D image
    tomo = load(fn);
    tomo.param_combine.array_param.freq_rng = new_freq_rng{fn_idx};
    
    % Load layer data
    param = tomo.param_combine;
    
    %% Get the generic layer data path
    layer_path = fullfile(ct_filename_out(param,param.combine.layer_fn,'',0));
    
    %% Load the current frame
    layer_fn = fullfile(layer_path,sprintf('Data_%s_%03d.mat',param.day_seg,param.combine.frm));
    layer = load(layer_fn);
    layer_gps_time = layer.GPS_time;
    surface = layer.layerData{1}.value{2}.data;
    bottom = layer.layerData{2}.value{2}.data;
    %% Get the previous frame if necessary
    if tomo.GPS_time(1) < layer_gps_time(1)-1
      layer_fn = fullfile(layer_path,sprintf('Data_%s_%03d.mat',param.day_seg,param.combine.frm-1));
      if exist(layer_fn,'file')
        layer = load(layer_fn);
        layer_gps_time = [layer.GPS_time layer_gps_time];
        surface = [layer.layerData{1}.value{2}.data surface];
        bottom = [layer.layerData{1}.value{2}.data bottom];
      end
    end
    %% Get the next frame if necessary
    if tomo.GPS_time(end) > layer_gps_time(end)+1
      layer_fn = fullfile(layer_path,sprintf('Data_%s_%03d.mat',param.day_seg,param.combine.frm+1));
      if exist(layer_fn,'file')
        layer = load(layer_fn);
        layer_gps_time = [layer_gps_time layer.GPS_time];
        surface = [surface layer.layerData{1}.value{2}.data];
        bottom = [bottom layer.layerData{1}.value{2}.data];
      end
    end
    %% Since layer files may have overlapping data, sort it
    [layer_gps_time new_surface_idxs] = sort(layer_gps_time);
    surface = surface(new_surface_idxs);
    bottom = bottom(new_surface_idxs);
    
    %% Do the interpolation
    good_idxs = find(isfinite(surface));
    tomo.param_combine.array_param.surface = interp1(layer_gps_time(good_idxs),surface(good_idxs),tomo.GPS_time,'linear','extrap');
    good_idxs = find(isfinite(bottom));
    tomo.param_combine.array_param.bottom = interp1(layer_gps_time(good_idxs),bottom(good_idxs),tomo.GPS_time,'linear','extrap');
    
    array_param = tomo.param_combine.array_param;
    
    tomo.valR ...
      = NaN*zeros(array_param.Nsv,size(tomo.Topography,3),'single');
    tomo.bins ...
      = NaN*zeros(array_param.Nsv,size(tomo.Topography,3),'single');
    tomo.val ...
      = NaN*zeros(size(tomo.Topography,1), size(tomo.Topography,3),'single');
    tomo.freq ...
      = NaN*zeros(size(tomo.Topography,1),size(tomo.Topography,3),'single');
    
    theta = array_param.sv_fh(array_param.Nsv,array_param.wfs.fc);
    theta_fftshift = fftshift(theta);
    
    array_param.freq_rng = array_param.freq_rng(array_param.freq_rng ~= 1);
    for rline = 1:size(tomo.Topography,3)
      
      % Determine the flat surface position
      surface = array_param.surface(rline) * c/2;
      thickness = (array_param.bottom(rline) - array_param.surface(rline)) * c/2/sqrt(3.15);
      
      bottom_range = surface./cos(theta_fftshift/1.3) ...
        + thickness./cos(asin(sin(theta_fftshift/0.5)/sqrt(3.15)));
      
      time_surf = array_param.surface(rline);
      range = zeros(size(array_param.wfs.time));
      range(array_param.wfs.time<time_surf) ...
        = array_param.wfs.time(array_param.wfs.time<time_surf) * c/2;
      range(array_param.wfs.time>time_surf) = surface ...
        + (array_param.wfs.time(array_param.wfs.time>time_surf)-time_surf) * (c/2/sqrt(er_ice));
      range_guard = [-100 100];
      
      bottom_near = round(interp1(range,1:length(range),bottom_range+range_guard(1)));
      
      bottom_far = round(interp1(range,1:length(range),bottom_range+range_guard(2)));
      bottom_far(isnan(bottom_far)) = length(range);
      
      if 0
        figure(1); clf;
        imagesc(lp(Sarray))
        hold on
        plot((1:64)-2,bottom_near,'k')
        plot((1:64)+2,bottom_far,'k')
        hold off
      end
      
      Sarray = tomo.Topography(:,:,rline);
      
      % Search for maximum value in each cross-track frequency bin
      % and keep track of which range bin the max falls in
      for freq_idx = 2:size(Sarray,2)
        [tomo.valR(freq_idx,rline) tomo.bins(freq_idx,rline)] = max(Sarray(bottom_near(freq_idx):bottom_far(freq_idx),freq_idx));
        tomo.bins(freq_idx,rline) = tomo.bins(freq_idx,rline) + bottom_near(freq_idx) - 1;
      end
      % Search for maximum value in each range bin and keep track of
      % which cross-track frequency the max falls in
      
      if swath_mode{fn_idx} == 1
        % Side looking
        %bottom_rbins = min(bottom_near(array_param.freq_rng)):max(bottom_far(array_param.freq_rng));
        bottom_rbins = min(bottom_near(array_param.freq_rng)) + (1:2000);
        bottom_freq_rng_near = round(interp1(bottom_near(array_param.freq_rng),array_param.freq_rng,bottom_rbins));
        bottom_freq_rng_far = round(interp1(bottom_far(array_param.freq_rng),array_param.freq_rng,bottom_rbins));
        min_freq = min(array_param.freq_rng);
        max_freq = max(array_param.freq_rng);
        if mean(array_param.freq_rng) < length(theta_fftshift)/2
          % Left looking
          bottom_freq_rng_near(isnan(bottom_freq_rng_near)) = min_freq;
          bottom_freq_rng_far(isnan(bottom_freq_rng_far)) = max_freq;
          for bottom_rbins_idx = 1:length(bottom_rbins)
            bin_idx = bottom_rbins(bottom_rbins_idx);
            [tomo.val(bin_idx,rline) tomo.freq(bin_idx,rline)] = max(Sarray(bin_idx,max(1,bottom_freq_rng_near(bottom_rbins_idx)-2):bottom_freq_rng_far(bottom_rbins_idx)+2));
            tomo.freq(bin_idx,rline) = tomo.freq(bin_idx,rline) + bottom_freq_rng_near(bottom_rbins_idx)-2 + 1;
          end
        else
          % Right looking
          bottom_freq_rng_near(isnan(bottom_freq_rng_near)) = max_freq;
          bottom_freq_rng_far(isnan(bottom_freq_rng_far)) = min_freq;
          for bottom_rbins_idx = 1:length(bottom_rbins)
            bin_idx = bottom_rbins(bottom_rbins_idx);
            [tomo.val(bin_idx,rline) tomo.freq(bin_idx,rline)] = max(Sarray(bin_idx,bottom_freq_rng_far(bottom_rbins_idx)-2:min(size(Sarray,2),bottom_freq_rng_near(bottom_rbins_idx)+2)));
            tomo.freq(bin_idx,rline) = tomo.freq(bin_idx,rline) + bottom_freq_rng_far(bottom_rbins_idx)-2 + 1;
          end
        end
      else
        % Nadir looking
        [tomo.val(:,rline) tomo.freq(:,rline)] = max(Sarray(:,array_param.freq_rng),[],2);
        tomo.freq(:,rline) = array_param.freq_rng(tomo.freq(:,rline));
      end
      
    end
    
    
    % Coefficients from tomography_interp.m
    Cy = 1.0e+03 * [-0.001660711272909
      0.090493698716517
      0.001048136114533
      0.000217252853224
      -0.541884104307672
      4.744743290064101
      0.003971414253168
      1.309756959638642
      -0.013945276404631
      0.699992786161553
      -0.032765372052905
      0.012051703640593
      -0.001483139287919
      0.424229484875972
      -0.002626573153978
      0.699992786161575
      -0.860590076181269
      -1.524460881137507
      -0.335553751187282
      0.588158785381890
      0.903908426900442
      0.053668026221286];
    Cz = 1.0e+05 * [3.508380657395748
      0.078304143910234
      3.510483542206905
      -0.011025811415163
      -0.001674155930760
      1.754974750327560
      -0.000289369503542
      0.000151185549539
      0.000429567889924
      -0.000548543415327
      0.590125520756759
      0.000283533742311
      -0.013579345401646
      -0.005777378623296
      0.002278003156437
      -3.510565644678335
      0.003606277246845
      0.065714846428169
      0.004452579897103
      0.065714846427901
      -0.007257256182412
      0.045070547933172];
    scale_factors = 1.0e+03 * [0.000000095286139   0.001047197551197   1.950000000000000];
    
    
    Surface = tomo.param_combine.array_param.surface;
    
    
    if swath_mode{fn_idx} == 1
      tomo.freq(lp(tomo.val) < 15) = NaN;
      %     freq_med = medfilt2(tomo.freq,[11 5]);
      %     bad_mask = abs(freq_med - tomo.freq) > 5;
      %     bad_mask = filter2(ones(11),bad_mask);
      %
      %     for rline = 1:size(tomo.Topography,3)
      %       tomo.freq(1:tomo.bins(33,rline),rline) = NaN;
      %       tomo.freq(tomo.bins(33,rline)+950:end,rline) = NaN;
      %     end
      %     tomo.freq(bad_mask > 0) = NaN;
      %
      % =====================================================================
      % Use polynomial to geocode results
      body_idx6_y = zeros(size(tomo.freq),'single');
      body_idx6_z = zeros(size(tomo.freq),'single');
      
      freq_filt = filter2(ones(5,3)/15, tomo.freq);
      
      for rbin = 1:size(freq_filt,1)
        timeDelay = array_param.wfs.time(rbin);
        
        theta_row = NaN*ones(size(freq_filt(rbin,:)));
        
        theta_row(~isnan(freq_filt(rbin,:))) = interp1(1:length(theta_fftshift), theta_fftshift, freq_filt(rbin,~isnan(freq_filt(rbin,:))));
        Y = abs(theta_row).'/scale_factors(2);
        X = timeDelay*ones(size(Y))/scale_factors(1);
        Z = (Surface.'*c/2)/scale_factors(3);
        
        body_idx6_y(rbin,:) = sign(theta_row.') .* ...
          ([ones(size(X)) Y Z X.^2 Y.^2 X.*Y X.*Z Y.*Z X.^3 Y.^3 X.^2.*Y X.^2.*Z Y.^2.*X Y.^2.*Z Z.^2.*X Y.^3 Y.^3.*X Y.^4 Y.^4.*Z Y.^5 Y.^5.*Z Y.^5.*X]*Cy);
        body_idx6_z(rbin,:) = ...
          [ones(size(X)) X Y Z X.^2 Y.^2 X.*Y X.*Z Y.*Z X.^3 Y.^3 X.^2.*Y Y.^2.*X Y.^2.*Z exp(X) exp(Y) X.*Y.^4 Y.^4 Z.*Y.^4 Y.^4 Y.^5.*Z Y.^5]*Cz;
      end
      
    else
      
      bad_mask = isnan(tomo.bins);
      tomo.bins(bad_mask) = 1;
      timeDelay = array_param.wfs.time(tomo.bins);
      timeDelay(bad_mask) = NaN;
      
      % =====================================================================
      % Use polynomial to geocode results
      body_idx6_y = zeros(size(timeDelay),'single');
      body_idx6_z = zeros(size(timeDelay),'single');
      for inc_idx = 1:size(timeDelay,1)
        X = timeDelay(inc_idx,:).'/scale_factors(1);
        Y = abs(theta_fftshift(inc_idx)*ones(size(X)))/scale_factors(2);
        Z = (Surface.'*c/2)/scale_factors(3);
        
        body_idx6_y(inc_idx,:) = sign(theta_fftshift(inc_idx)) .* ...
          [ones(size(X)) Y Z X.^2 Y.^2 X.*Y X.*Z Y.*Z X.^3 Y.^3 X.^2.*Y X.^2.*Z Y.^2.*X Y.^2.*Z Z.^2.*X Y.^3 Y.^3.*X Y.^4 Y.^4.*Z Y.^5 Y.^5.*Z Y.^5.*X]*Cy;
        body_idx6_z(inc_idx,:) = ...
          [ones(size(X)) X Y Z X.^2 Y.^2 X.*Y X.*Z Y.*Z X.^3 Y.^3 X.^2.*Y Y.^2.*X Y.^2.*Z exp(X) exp(Y) X.*Y.^4 Y.^4 Z.*Y.^4 Y.^4 Y.^5.*Z Y.^5]*Cz;
      end
      
    end
    
    % [B,A] = butter(2,0.01);
    % tomo.param_combine.array_param.fcs{1}{1}.y(1,:) = filtfilt(B,A,tomo.param_combine.array_param.fcs{1}{1}.y(1,:));
    % tomo.param_combine.array_param.fcs{1}{1}.y(2,:) = filtfilt(B,A,tomo.param_combine.array_param.fcs{1}{1}.y(2,:));
    % tomo.param_combine.array_param.fcs{1}{1}.y(3,:) = filtfilt(B,A,tomo.param_combine.array_param.fcs{1}{1}.y(3,:));
    % tomo.param_combine.array_param.fcs{1}{1}.z(1,:) = filtfilt(B,A,tomo.param_combine.array_param.fcs{1}{1}.z(1,:));
    % tomo.param_combine.array_param.fcs{1}{1}.z(2,:) = filtfilt(B,A,tomo.param_combine.array_param.fcs{1}{1}.z(2,:));
    % tomo.param_combine.array_param.fcs{1}{1}.z(3,:) = filtfilt(B,A,tomo.param_combine.array_param.fcs{1}{1}.z(3,:));
    
    
    x_plane = NaN*zeros(size(body_idx6_y));
    y_plane = NaN*zeros(size(body_idx6_y));
    z_plane = NaN*zeros(size(body_idx6_y));
    
    body_idx6_z = body_idx6_z + repmat(Surface,[size(body_idx6_z,1) 1]) * c/2;
    
    for rline = 1:size(body_idx6_y,2)
      x_plane(:,rline) = tomo.param_combine.array_param.fcs{1}{1}.origin(1,rline) ...
        + tomo.param_combine.array_param.fcs{1}{1}.y(1,rline) * body_idx6_y(:,rline) ...
        + tomo.param_combine.array_param.fcs{1}{1}.z(1,rline) * -body_idx6_z(:,rline);
      y_plane(:,rline) = tomo.param_combine.array_param.fcs{1}{1}.origin(2,rline) ...
        + tomo.param_combine.array_param.fcs{1}{1}.y(2,rline) * body_idx6_y(:,rline) ...
        + tomo.param_combine.array_param.fcs{1}{1}.z(2,rline) * -body_idx6_z(:,rline);
      z_plane(:,rline) = tomo.param_combine.array_param.fcs{1}{1}.origin(3,rline) ...
        + tomo.param_combine.array_param.fcs{1}{1}.y(3,rline) * body_idx6_y(:,rline) ...
        + tomo.param_combine.array_param.fcs{1}{1}.z(3,rline) * -body_idx6_z(:,rline);
    end
    
    [points.lat,points.lon,points.elev] = ecef2geodetic(x_plane,y_plane,z_plane,WGS84.ellipsoid);
    points.lat = points.lat * 180/pi;
    points.lon = points.lon * 180/pi;
    points.swath_mode = swath_mode{fn_idx};
    points.freq_rng = new_freq_rng{fn_idx};
    points.y = body_idx6_y;
    
    if ~exist(base_path2,'dir')
      mkdir(base_path2)
    end
    [~,fn_name] = fileparts(fn);
    out_fn = fullfile(base_path2,fn_name);
    fprintf('  Saving %s\n', out_fn);
    save(out_fn,'points');
    
    figure(1); clf;
    imagesc(points.elev);
    drawnow;
    pause(1);
    
  end
  return
  
  
  
  
  
  
  
  
  
  
  
  
  
  
else
  
  
  tic;
  physical_constants;
  
  % =====================================================================
  % User Settings
  % =====================================================================
  
  fns = {};
  new_freq_rng = {};
  swath_mode = {};
  
  base_path2 = '/mnt/products/rds/2013_Antarctica_Basler/CSARP_tomo/20131220_03';
  for frm = [3 4 5 6]
    fns{end+1} = fullfile(base_path2,sprintf('Data_img_01_20131220_03_%03i',frm));
    new_freq_rng{end+1} = [1:64];
    swath_mode{end+1} = 2;
    fns{end+1} = fullfile(base_path2,sprintf('Data_img_02_20131220_03_%03i',frm));
    new_freq_rng{end+1} = [33:64];
    swath_mode{end+1} = 1;
    fns{end+1} = fullfile(base_path2,sprintf('Data_img_03_20131220_03_%03i',frm));
    new_freq_rng{end+1} = [1:33];
    swath_mode{end+1} = 1;
  end
  
  x = [];
  y = [];
  z = [];
  lat = [];
  lon = [];
  y_offset = [];
  swath_idx = [];
  swath_idx_maps = [];
  for fn_idx = 1:length(fns)
    fn = fns{fn_idx};
    
    fprintf('Reading %s\n', fn);
    
    %% Load 3D image
    load(fn);
    
    %% Determine the swath for this file
    if swath_mode{fn_idx} == 1
      if new_freq_rng{fn_idx}(1) == 1
        good_mask = points.y < -250 & points.y > -1250;
      else
        good_mask = points.y > 250 & points.y < 1250;
      end
    else
      good_mask = points.y > -550 & points.y < 550;
    end
    
%     ff = points.elev;
%     ff(~good_mask) = NaN;
%     figure(1); clf;
%     imagesc(ff);
%     caxis([-700 -450])
%     keyboard
    [tmp_x,tmp_y] = geodetic_to_stereographic(points.lat(good_mask),points.lon(good_mask));
    x = [x reshape(tmp_x,[1 numel(tmp_x)])];
    y = [y reshape(tmp_y,[1 numel(tmp_y)])];
    lat = [lat reshape(points.lat(good_mask),[1 numel(points.lat(good_mask))])];
    lon = [lon reshape(points.lon(good_mask),[1 numel(points.lon(good_mask))])];
    z = [z reshape(points.elev(good_mask),[1 numel(points.elev(good_mask))])];
    y_offset = [y_offset reshape(points.y(good_mask),[1 numel(points.y(good_mask))])];
    swath_idx = [swath_idx fn_idx*ones([1 numel(points.y(good_mask))])];
  end
  
  y_offset_min = min(y_offset);
  y_offset_max = max(y_offset);
  p = polyfit(y_offset,z,4);
  
  figure(1); clf;
  plot(y_offset,z,'.');
  hold on;
  plot(y_offset_min:y_offset_max, polyval(p,y_offset_min:y_offset_max),'r');
  hold off;
  %
  z = z - (polyval(p,y_offset) - polyval(p,0));
  %
  % figure(2); clf;
  % plot(y_offset,z,'.');
  %
  % return
  
  %
  % good_mask = ~isnan(x) & ~isnan(y) & ~isnan(z);
  % z_grid = griddata(x(good_mask),y(good_mask),z(good_mask),x_axis,y_axis.');
  
  good_mask = ~isnan(x) & ~isnan(y) & ~isnan(z);
  x = x(good_mask);
  y = y(good_mask);
  z = z(good_mask);
  
  x_min = min(x)-1;
  x_max = max(x)+1;
  y_min = min(y)-1;
  y_max = max(y)+1;
  
  dx = 25; dy = 25;
  x_axis = x_min:dx:x_max;
  y_axis = y_min:dy:y_max;
  
  z_grid = NaN*zeros(length(y_axis),length(x_axis));
  error_map = NaN*zeros(length(y_axis),length(x_axis));
  swath_idx_maps = zeros(length(y_axis),length(x_axis),length(fns),'uint8');
  
  dist_thresh = 50;
  % Determine the cell each point will be in.
  cell_x = floor((x-x_min) / dx);
  cell_y = floor((y-y_min) / dy);
  
  for x_idx = 1:length(x_axis)
    fprintf('%d of %d (%s)\n', x_idx, length(x_axis),datestr(now));
    x_sub = x(cell_x == x_idx-1);
    y_sub = y(cell_x == x_idx-1);
    z_sub = z(cell_x == x_idx-1);
    cell_x_sub = cell_x(cell_x == x_idx-1);
    cell_y_sub = cell_y(cell_x == x_idx-1);
    for y_idx = 1:length(y_axis)
      z_grid(y_idx,x_idx) = mean(z_sub(cell_x_sub == x_idx-1 & cell_y_sub == y_idx-1));
      error_map(y_idx,x_idx) = std(z_sub(cell_x_sub == x_idx-1 & cell_y_sub == y_idx-1));
      swath_idx_maps(y_idx,x_idx,swath_idx(cell_x_sub == x_idx-1 & cell_y_sub == y_idx-1)) = 1;
      %       dist = sqrt(abs(x_axis(x_idx) - x).^2 + abs(y_axis(y_idx) - y).^2);
      %       good_idxs = dist < dist_thresh;
      %       if any(good_idxs)
      %         weights = dist_thresh - dist(good_idxs);
      %         z_grid(y_idx,x_idx) = sum(z(good_idxs).*weights/sum(weights));
      %       end
    end
  end
  
  
  save('/mnt/products/test_tomography2.mat','z_grid','x_axis','y_axis','error_map','swath_idx_maps');
  
  return
end

save('/mnt/products/test_tomography3.mat','z_grid','x_axis','y_axis','error_map','swath_idx_maps');
