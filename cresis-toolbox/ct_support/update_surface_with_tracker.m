% script update_surface_with_tracker
%
% Updates the Surface variable in an echogram using one of three surface
% detections methods.

%% User Settings
% ----------------------------------------------------------------------
% params = read_param_xls(ct_filename_param('snow_param_2009_Greenland_P3.xls'),'20090421_04','post');
params = read_param_xls(ct_filename_param('snow_param_2010_Greenland_DC8.xls'),'20100419_02','post');
% params = read_param_xls(ct_filename_param('snow_param_2011_Greenland_P3.xls'),'20110325_02','post');
% params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'),'','post');
params.cmd.generic = 1;
% params.cmd.frms = [140:520];
% params.cmd.frms = [1:5]; %0412_02
% params.cmd.frms = [1:113];%0412_3
params.cmd.frms = [131:360]; %419_02
% params.cmd.frms = [228:521]; %20_02
% params.cmd.frms = [215:229,383:394]; %21_01
% params.cmd.frms = [24:28,128:132]; %21_03
debug_level = 0;


echogram_source = 'deconv'; % location of echogram data
layerdata_source = 'layerData'; % location of layerData
% save_sources: cell array of strings indicating which layer sources
%   should be updated (options are ops, layerData, and echogram)
save_sources = {'echogram'};
save_ops_layer = 'surface';

surf_override = [];
surf_override.en = true;
surf_override.min_bin = 0;
surf_override.manual = false;

% 
% % RDS threshold
% surf_override.method = 'threshold';
% surf_override.noise_rng = [0 -50 -10];
% surf_override.min_bin = 0.75e-6;
% surf_override.threshold = 9;
% surf_override.sidelobe	= 15;
% surf_override.max_diff	= inf;
% surf_override.filter_len	= 1;
% surf_override.medfilt = 3;
% surf_override.search_rng	= 0:2;

% % RDS
% surf_override.method = 'max';
% surf_override.min_bin = 0.75e-6;
% surf_override.search_rng	= [-11:0];
% surf_override.threshold = 13;
% surf_override.medfilt = 3;
% FEEDTHRU IS SEASON DEPENDENT (this example is for RDS 2009 Greenland TO)
% surf_override.feedthru.time = 1e-9*[102;1069;1803;2436;3003;3603];
% surf_override.feedthru.power_dB = [-50.9;-82.0;-106.3;-128.9;-131.7;-152.1];

% % Accum
% surf_override.method = 'threshold';
% surf_override.noise_rng = [200 -300 -100];
% surf_override.min_bin = 0.1e-6;
% surf_override.threshold = 9;
% surf_override.sidelobe	= 12;
% surf_override.max_diff	= inf;
% surf_override.filter_len	= 5;
% surf_override.search_rng	= 0:1;

% % FMCW Land Ice
surf_override.method = 'threshold';
surf_override.noise_rng = [100 -400 -200];
surf_override.min_bin = 1e-6;
surf_override.threshold = 7;
surf_override.sidelobe	= 13;
surf_override.max_diff	= 30e-9;
surf_override.filter_len	= [3 13];
surf_override.search_rng	= 0:30;
surf_override.detrend = 0;
% surf_override.init.method	= 'medfilt';
% surf_override.init.medfilt	= 11;
surf_override.init.method	= 'dem';
surf_override.init.dem_offset = 30e-9;
% surf_override.init.method	= 'snake';
% surf_override.init.search_rng	= [-240:240];
% surf_override.method = 'snake';
% surf_override.search_rng	= -90:90;

% 
% % FMCW Sea Ice
% surf_override.method = 'threshold';
% surf_override.noise_rng = [100 -400 -100];
% surf_override.threshold = 9;
% surf_override.sidelobe	= 13;
% surf_override.max_diff	= 30e-9;
% surf_override.filter_len	= 7;
% surf_override.search_rng	= 0:30;
% surf_override.detrend = 2;
% surf_override.init.method	= 'medfilt';
% surf_override.init.medfilt	= 51;

%% Automated Section
% ----------------------------------------------------------------------
physical_constants;

% gimp_geoid_loaded = false;

for param_idx = 1:length(params)
  param = params(param_idx);
  
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  orig_surf = merge_structs(param.get_heights.surf,surf_override);
  
  if ~orig_surf.en
    continue;
  end
  
  fprintf('Updating surface %s (%s)\n', param.day_seg, datestr(now,'HH:MM:SS'));
  
  %% Load in GIMP and Geoid
  if strcmpi(orig_surf.init.method,'dem') & ~gimp_geoid_loaded
    % Load in Geoid
    fn = ct_filename_gis([],'world\egm96_geoid\WW15MGH.DAC');
    [egm96_lat,egm96_lon,egm96] = egm96_loader(fn);
    egm96_lon = [egm96_lon 360];
    egm96 = [egm96 egm96(:,1)];
    
    % Interpolate to find GIMP elevation along flight path
    [gimp_dem, gimp_R, tmp] = geotiffread('/cresis/snfs1/dataproducts/GIS_data/greenland/DEM/GIMP/gimpdem_90m.tif');
    gimp_dem = double(gimp_dem);
    gimp_proj = geotiffinfo('/cresis/snfs1/dataproducts/GIS_data/greenland/DEM/GIMP/gimpdem_90m.tif');
    gimp_geoid_loaded = true;
  end
  
  %% Load "frames" for this segment
  load(ct_filename_support(param,'','frames'));
  
  %% Determine valid frames to process
  if isempty(param.cmd.frms)
    param.cmd.frms = 1:length(frames.frame_idxs);
  end
  % Remove frames that do not exist from param.cmd.frms list
  [valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
  if length(valid_frms) ~= length(param.cmd.frms)
    bad_mask = ones(size(param.cmd.frms));
    bad_mask(keep_idxs) = 0;
    warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
      param.cmd.frms(find(bad_mask,1)));
    param.cmd.frms = valid_frms;
  end

  %% Track each of the frames
  for frm_idx = 1:length(param.cmd.frms)
    frm = param.cmd.frms(frm_idx);
    
    data_fn = fullfile(ct_filename_out(param,echogram_source,''), ...
      sprintf('Data_%s_%03d.mat', param.day_seg, frm));
    fprintf('  Updating %s (%s)\n', data_fn, datestr(now,'HH:MM:SS'));
    
    if ~exist(data_fn)
      fprintf('    Missing file\n');
      continue;
    end
    
    %% Load echogram data
    mdata = load(data_fn);
    
    %% Update surf structure for this particular echogram
    %   - min_bin field is specified as two way travel time, but needs to
    %   be converted to range bins
    surf = orig_surf;
    surf.min_bin = find(mdata.Time > orig_surf.min_bin, 1);
    if isfield(orig_surf,'max_bin')
      surf.max_bin = find(mdata.Time > orig_surf.max_bin, 1);
    else
      surf.max_bin = inf;
    end
    dt = mdata.Time(2) - mdata.Time(1);

    if ~isfield(surf,'manual')
      surf.manual = false;
    end
    
    if isfield(surf,'feedthru')
      %% Optional feed through removal
      
      % Interpolate feed through power levels on to data time axis
      feedthru_threshold = interp1(surf.feedthru.time,surf.feedthru.power_dB,mdata.Time);
      feedthru_threshold = interp_finite(feedthru_threshold);
      
      % Set all data to zero that does not exceed the feed through
      % threshold power
      for rline=1:size(mdata.Data,2)
        mdata.Data(:,rline) = mdata.Data(:,rline) .* (lp(mdata.Data(:,rline)) > feedthru_threshold);
      end
    end
    
    %% Interpolate GIMP and Geoid
    if strcmpi(orig_surf.init.method,'dem')
      mdata.geoid_elev = interp2(egm96_lon,egm96_lat,egm96,mod(mdata.Longitude,360),mdata.Latitude);
      [mdata.x,mdata.y] = projfwd(gimp_proj,mdata.Latitude,mdata.Longitude);
      mdata.gimp_dem = interp2(gimp_dem,(mdata.x-gimp_R(3,1))/gimp_R(2,1)+1,(mdata.y-gimp_R(3,2))/gimp_R(1,2)+1);
      
      % Merge GIMP and Geoid
      surf.dem = mdata.gimp_dem;
      surf.dem(isnan(surf.dem)) = mdata.geoid_elev(isnan(surf.dem));
      surf.dem = (mdata.Elevation - surf.dem) / (c/2);
      surf.dem = interp1(mdata.Time,1:length(mdata.Time),surf.dem + surf.init.dem_offset);
      surf.dem = interp_finite(surf.dem,length(mdata.Time)/2);
    end
    
    surf.max_diff = orig_surf.max_diff/dt;
    
    %% Track the surface
    if surf.manual
      [new_surface,pnt] = tracker_snake_simple(mdata.Data,surf);
      fprintf('  Press F1 for help\n');
      layer = tracker_snake_manual_gui(lp(mdata.Data),pnt);

    elseif strcmpi(surf.method,'threshold')
      new_surface = tracker_threshold(mdata.Data,surf);
    elseif strcmpi(surf.method,'max')
      new_surface = tracker_max(mdata.Data,surf);
    elseif strcmpi(surf.method,'snake')
      new_surface = tracker_snake_simple(mdata.Data,surf);
    else
      error('Not a supported surface tracking method.');
    end

    %% Run median filtering on tracked surface
    if isfield(surf,'medfilt') && ~isempty(surf.medfilt)
      % OLD METHOD: new_surface = medfilt1(new_surface,surf.medfilt);
      surf.medfilt_threshold = 10;
      for rline=1:length(new_surface)
        rlines = rline + (-surf.medfilt:surf.medfilt);
        rlines = rlines(rlines>=1 & rlines<=length(new_surface));
        if abs(new_surface(rline) - nanmedian(new_surface(rlines))) > surf.medfilt_threshold
          new_surface(rline) = nanmedian(new_surface(rlines));
        end
      end
    end
    
    % Convert from range bins to two way travel time
    Surface = interp1(1:length(mdata.Time), mdata.Time, new_surface);

    if debug_level > 0
      % Debug result
      figure(1); clf;
      imagesc([],mdata.Time,lp(mdata.Data));
      colormap(1-gray(256));
      hold on;
      plot(Surface);
      if strcmpi(orig_surf.init.method,'dem')
        plot(interp1(1:length(mdata.Time),mdata.Time,surf.dem),'g')
        plot(interp1(1:length(mdata.Time),mdata.Time,surf.dem-surf.max_diff),'r')
        plot(interp1(1:length(mdata.Time),mdata.Time,surf.dem+surf.max_diff),'b')
      end
      hold off;
      ylim([min(Surface)-10e-9 max(Surface)+10e-9])
      keyboard
    end
    
    %% Save result
    if any(strcmpi(save_sources,'echogram'))
      save(data_fn,'-append','Surface');
    end
    
    if any(strcmpi(save_sources,'layerdata'))
      layer_fn = fullfile(ct_filename_out(param,layerdata_source,''), ...
        sprintf('Data_%s_%03d.mat', param.day_seg, frm));
      if exist(layer_fn,'file')
        % Load the layerData file
        lay = load(layer_fn);
        % Update the surface auto picks
        lay.layerData{1}.value{2}.data = interp1(mdata.GPS_time,Surface,lay.GPS_time);
        lay.layerData{1}.value{2}.data = interp_finite(lay.layerData{1}.value{2}.data);
        % Append the new results back to the layerData file
        save(layer_fn,'-append','-struct','lay','layerData');
      end
    end
    
    if any(strcmpi(save_sources,'ops'))
      % OPS query to get the point path ID's
      ops_param = struct('properties',[]);
      ops_param.properties.location = param.post.ops.location;
      ops_param.properties.season = param.season_name;
      ops_param.properties.segment = param.day_seg;
      ops_param.properties.start_gps_time = mdata.GPS_time(1);
      ops_param.properties.stop_gps_time = mdata.GPS_time(end);
      ops_param.properties.lyr_name = save_ops_layer;
      
      sys = ct_output_dir(param.radar_name);
      [status,data] = opsGetLayerPoints(sys,ops_param);
      
      % Write the new layer information to these point path ID's
      ops_param = struct('properties',[]);
      ops_param.properties.point_path_id = data.properties.point_path_id;
      ops_param.properties.twtt = interp1(mdata.GPS_time,Surface,data.properties.gps_time);
      ops_param.properties.type = 2*ones(size(ops_param.properties.twtt));
      ops_param.properties.quality = 1*ones(size(ops_param.properties.twtt));
      ops_param.properties.lyr_name = save_ops_layer;
      
      opsCreateLayerPoints(sys,ops_param);
    end
  end
  
end






