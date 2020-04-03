function create_surfdata(param,mdata)
% tomo.create_surfdata(param,mdata)
%
% Description: Usually this function is called from tomo_collate_task.
%   Using a surface DEM and an ice mask, this function adds an aligned
%   surface dem and ice mask to a file.
%
% Inputs:
% param: struct from parameter spreadsheet
%  .tomo_collate: struct which control create_surfdata
%   .out_dir: ct_filename_out directory to which output
%     surfData will be exported. Default is surfData.
%   .layer_params: opsLoadLayers parameter struct array. The first element
%     should load the ice top (surface) layer and the second element
%     should load the ice bottom.
%   .surfData_mode: surfData mode string, can be one of:
%     'overwrite': will overwrite the whole file
%     'append': will only update the specified layers
%     'fillgaps': will only update the specified layers if they did not
%       already exist
%   .surfdata_cmds: struct of surfdata commands to run
%     .cmd: String with command: 'detect', 'dem', 'extract', 'viterbi', or 'trws'
%     .surf_names: String or cell array of strings with surface names to be
%       updated by the results of the command
%     .visible: The visibility setting for this layer
%     .plot_name_values: Plot name-value property fields
%     DETECT parameters
%     .data_threshold: pixel values above this will be clipped (default is 13.5)
%     DEM parameters
%     .dem_fn: filename to geotiff with DEM
%     EXTRACT parameters
%     .data_threshold: pixel values above this will be clipped (default is 13.5)
%     VITERBI parameters
%     .smooth_weight: smoothness weight (scalar, default is 55)
%     .smooth_var: smoothness variance (scalar, default is -1)
%     .repulsion: surface repulsion scaling factor (scalar, default is 150)
%     .egt_weight: extra ground-truth weight (scalar, default is 10)
%     .ice_bin_thr: ice_mask scanning threshold (scalar, default is 3)
%     TRWS parameters
%     .smooth_weight: smoothness weight [] (default is [22 22])
%     .smooth_var: Gaussian weighting in the elevation angle bin dimension (default is 32)
%     .max_loops: number of loops to run for (default is 50)
%   .bounds_relative: DOA bins and along-track slices to trim off from each edge [top bottom left right]
%     E.g. [0 0 0 0] would use all pixels, [1 0 0 0] would ignore the first row.
%     Usually trimming a few rows off the top/bottom is a good idea if
%     these represent angles out to +/- 90 deg.
%
% mdata: 3D data file struct (Data_YYYYMMDD_SS_FFF.mat with Tomo)
%   field)
%
% Outputs:
%   NONE
%
% See also: tomo.run_collate, tomo.collate, tomo_collate_task,
%   tomo.fuse_images, tomo.add_icemask_surfacedem, tomo.create_surfdata,
%
% Author: John Paden, Jordan Sprick, Mingze Xu, and Victor Berger

if ~isfield(param.tomo_collate,'surf_out_path') || isempty(param.tomo_collate.surf_out_path)
  param.tomo_collate.surf_out_path = 'surfData';
end

% If DOA method is used, set doa_method_flag = true
array_proc_methods; % This script assigns the integer values for each method
if ischar(param.array.method)
  % Convert array method string to integer
  method_integer = [];
  if regexpi(param.array.method,'music_doa')
    method_integer(end+1) = MUSIC_DOA_METHOD;
  end
  if regexpi(param.array.method,'mle')
    method_integer(end+1) = MLE_METHOD;
  end
  if regexpi(param.array.method,'dcm')
    method_integer(end+1) = DCM_METHOD;
  end
  %   if regexpi(param.array.method,'pf')
  %     method_integer(end+1) = PF_METHOD;
  %   end
end
method_integer = intersect(method_integer, ...
  [MUSIC_DOA_METHOD MLE_METHOD DCM_METHOD PF_METHOD], 'stable');
if ~isempty(method_integer)
  doa_method_flag = true;
else
  doa_method_flag = false;
end

if ~isfield(param.tomo_collate,'merge_bottom_above_top') ...
    || isempty(param.tomo_collate.merge_bottom_above_top)
  param.tomo_collate.merge_bottom_above_top = 1;
end
merge_bottom_above_top = param.tomo_collate.merge_bottom_above_top;

%% Load surface and bottom information
param_load_layers = param;
param_load_layers.cmd.frms = round([-1,0,1] + param.load.frm);

layers = opsLoadLayers(param_load_layers,param.tomo_collate.layer_params);

%% Interpolate surface and bottom information to mdata
master = [];
master.GPS_time = mdata.GPS_time;
master.Latitude = mdata.Latitude;
master.Longitude = mdata.Longitude;
master.Elevation = mdata.Elevation;
for lay_idx = 1:length(layers)
  ops_layer = [];
  ops_layer{1}.gps_time = layers(lay_idx).gps_time;
  
  ops_layer{1}.type = layers(lay_idx).type;
  ops_layer{1}.quality = layers(lay_idx).quality;
  ops_layer{1}.twtt = layers(lay_idx).twtt;
  ops_layer{1}.type(isnan(ops_layer{1}.type)) = 2;
  ops_layer{1}.quality(isnan(ops_layer{1}.quality)) = 1;
  lay = opsInterpLayersToMasterGPSTime(master,ops_layer,[300 60]);
  layers(lay_idx).twtt_ref = lay.layerData{1}.value{2}.data;
end
Surface = layers(1).twtt_ref;
Bottom = layers(2).twtt_ref;

if length(Bottom)~=size(mdata.Tomo.img,3)
  error('This should not happen.');
  Bottom = mdata.Bottom;
end
if length(Surface)~=size(mdata.Tomo.img,3)
  error('This should not happen.');
  Surface = mdata.Surface;
end

%% Interpolate Bottom, mdata.twtt from twtt to bins
if ~doa_method_flag
  Bottom_bin = interp1(mdata.Time, 1:length(mdata.Time), Bottom);
  Bottom_bin(isnan(Bottom_bin)) = -1;
end
if ~isfield(mdata,'twtt')
  mdata.twtt = layers(1).twtt;
end
if isfield(mdata,'ice_mask')
  ice_mask = mdata.ice_mask;
else
  ice_mask = ones(size(mdata.twtt));
end

%% Surface tracking prep: Convert img to double and log-scale
if ~doa_method_flag
  data = 10*log10(double(mdata.Tomo.img));
end

%% Surface tracking prep
% 1. Convert from twtt to bins
twtt_bin = round(interp1(mdata.Time, 1:length(mdata.Time), mdata.twtt));
% 2. The tracking software assumes that the surface never approaches
%    within mu_length of the top/bottom of the range line, so we truncate
%    surface to ensure this never happens.
mu_length = 11;
twtt_bin(isnan(twtt_bin) | twtt_bin > length(mdata.Time)-mu_length) = length(mdata.Time)-mu_length;
if doa_method_flag
  twtt_bin(isnan(mdata.twtt) | (mdata.twtt==0)) = NaN;
end

%% Create output filename
out_dir = ct_filename_out(param,param.tomo_collate.surf_out_path,'');
if ~isdir(out_dir)
  mkdir(out_dir);
end
out_fn_name = sprintf('Data_%s_%03.0f.mat',param.day_seg,param.load.frm);
out_fn = fullfile(out_dir,out_fn_name);
fprintf('  Updating %s\n', out_fn);
if any(strcmpi(param.tomo_collate.surfData_mode,{'append','fillgaps'}))
  if ~exist(out_fn,'file')
    warning('%s mode selected, but output file does not exist. Switching to "overwrite" mode.',param.tomo_collate.surfData_mode);
    param.tomo_collate.surfData_mode = 'overwrite';
  else
    try
      sd = tomo.surfdata(out_fn);
    catch ME
      % Output file is not good, so we need to create
      warning('Output surfData file exists, but could not be loaded. Run "dbcont" to overwrite the file.');
      keyboard
      param.tomo_collate.surfData_mode = 'overwrite';
    end
  end
elseif ~strcmpi(param.tomo_collate.surfData_mode,'overwrite')
  error('Invalid surfData_mode %s.', param.tomo_collate.surfData_mode);
end

if strcmpi(param.tomo_collate.surfData_mode,'overwrite')
  %% Create surfData
  sd = tomo.surfdata();
  sd.radar_name = mdata.param_array.radar_name;
  sd.season_name = mdata.param_array.season_name;
  sd.day_seg = mdata.param_array.day_seg;
  sd.frm = mdata.param_array.load.frm;
  sd.gps_time = mdata.GPS_time;
  if ~doa_method_flag
    sd.theta = mdata.Tomo.theta(:,1);
  else
    sd.theta = mdata.Tomo.theta;
  end
  sd.time = mdata.Time(:); % Make a column vector
  sd.FCS.origin = mdata.param_array.array_proc.fcs.origin;
  sd.FCS.x = mdata.param_array.array_proc.fcs.x;
  sd.FCS.y = mdata.param_array.array_proc.fcs.y;
  sd.FCS.z = mdata.param_array.array_proc.fcs.z;
end
      
if ~doa_method_flag
  Nsv = size(mdata.Tomo.img,2);
else
  Nx = size(mdata.Tomo.theta,3);
  for Nx_idx = 1:Nx
    theta_tmp = mdata.Tomo.theta(:,:,Nx_idx);
    theta_tmp = theta_tmp(~isnan(theta_tmp));
    max_Nsv(Nx_idx) = length(theta_tmp);
  end
  Nsv = max(max_Nsv);
end

if ~doa_method_flag
  % Beamforming method
  try
    surf = sd.get_surf('top');
    if strcmpi(param.tomo_collate.surfData_mode,'overwrite')
      surf.y = twtt_bin;
      sd.set_surf(surf);
    end
  catch ME
    surf = tomo.surfdata.empty_surf();
    surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
    surf.y = twtt_bin;
    surf.plot_name_values = {'color','black','marker','x'};
    surf.name = 'top';
    sd.insert_surf(surf);
  end
else
  % DOA method
  try
    surf = sd.get_surf('top');
    if strcmpi(param.tomo_collate.surfData_mode,'overwrite')
      surf.y = twtt_bin;
      sd.set_surf(surf);
    end
  catch ME
    surf = tomo.surfdata.empty_surf();
    surf.y = twtt_bin;
    surf.x = NaN(size(twtt_bin));
    for rline_idx = 1:Nx
      theta_rline = mdata.Tomo.theta(:,:,rline_idx);
      if ~all(isnan(theta_rline(:)))
        theta_rline = theta_rline(~isnan(theta_rline));
        surf.x(1:length(theta_rline),rline_idx) = theta_rline;
      end
    end
    % Sort DOA min to max (and, accordingly, range-bins). But surf.y is
    % already sorted inside add_icemask_surfacedem
    [surf.x x_idx] = sort(surf.x*180/pi,1,'ascend');
%     for rline_idx = 1:Nx
%       surf.y(:,rline_idx) = surf.y(x_idx(:,rline_idx),rline_idx);
%     end
    
    surf.plot_name_values = {'color','black','marker','*'}; % 'x'
    surf.name = 'top';
    sd.insert_surf(surf);
  end
  ice_top.x = surf.x;
  ice_top.y = surf.y;
end


if ~doa_method_flag
  % Beamforming method
  try
    surf = sd.get_surf('bottom');
    if strcmpi(param.tomo_collate.surfData_mode,'overwrite')
      surf.y = NaN(size(twtt_bin));
      sd.set_surf(surf);
    end
  catch ME
    surf = tomo.surfdata.empty_surf();
    surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
    surf.y = NaN(size(twtt_bin));
    surf.plot_name_values = {'color','blue','marker','^'}; 
    surf.name = 'bottom';
    sd.insert_surf(surf);
  end
else
  % DOA method
  try
    surf = sd.get_surf('bottom');
    if strcmpi(param.tomo_collate.surfData_mode,'overwrite')
      surf.y = NaN(size(twtt_bin));
      sd.set_surf(surf);
    end
  catch ME
    surf = tomo.surfdata.empty_surf();
    surf.x = NaN(size(twtt_bin));
    for rline_idx = 1:Nx
      theta_rline = mdata.Tomo.theta(:,:,rline_idx);
      theta_rline = theta_rline(~isnan(theta_rline));
      if ~all(isnan(theta_rline(:)))
        surf.x(1:length(theta_rline),rline_idx) = theta_rline;
      end
    end
    surf.y = NaN(size(twtt_bin)); % Will be created later in this  script
    % Sort DOA min to max (and, accordingly, range-bins)
    [surf.x x_idx] = sort(surf.x*180/pi,1,'ascend');
    for rline_idx = 1:Nx
      surf.y(:,rline_idx) = surf.y(x_idx(:,rline_idx),rline_idx);
    end
  % Ensure non-negative ice thickness
  if merge_bottom_above_top && exist('ice_top','var') && isfield(ice_top,'y') && ~isempty(ice_top.y)
    surf.y(surf.y<ice_top.y) = ice_top.y(surf.y<ice_top.y);
  end
    
    surf.plot_name_values = {'color','blue','marker','o'}; % '^'
    surf.name = 'bottom';
    sd.insert_surf(surf);
  end
  plot_name_values = surf.plot_name_values;
end

if ~doa_method_flag
  % Beamforming method
  try
    surf = sd.get_surf('ice mask');
    if strcmpi(param.tomo_collate.surfData_mode,'overwrite')
      surf.y = mdata.ice_mask;
      sd.set_surf(surf);
    end
  catch ME
    surf = tomo.surfdata.empty_surf();
    surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
    surf.y = mdata.ice_mask;
    surf.plot_name_values = {'color','white','marker','x'};
    surf.name = 'ice mask';
    sd.insert_surf(surf);
  end
else
  % DOA method
  try
    surf = sd.get_surf('ice mask');
    if strcmpi(param.tomo_collate.surfData_mode,'overwrite')
      surf.y = mdata.ice_mask;
      sd.set_surf(surf);
    end
  catch ME
    surf = tomo.surfdata.empty_surf();
    surf.x = NaN(size(twtt_bin));
    for rline_idx = 1:Nx
      theta_rline = mdata.Tomo.theta(:,:,rline_idx);
      if ~all(isnan(theta_rline(:)))
        theta_rline = theta_rline(~isnan(theta_rline));
        surf.x(1:length(theta_rline),rline_idx) = theta_rline;
      end
    end
        surf.y = mdata.ice_mask;
        % Sort DOA min to max (and, accordingly, range-bins). But surf.y is
        % already sorted inside add_icemask_surfacedem
        [surf.x x_idx] = sort(surf.x*180/pi,1,'ascend');
%     for rline_idx = 1:Nx
%       surf.y(:,rline_idx) = surf.y(x_idx(:,rline_idx),rline_idx);
%     end
%     surf.plot_name_values = {'color','white','marker','x'};
    surf.plot_name_values = {'color',[0 0 0.5],'marker','x'};
    surf.name = 'ice mask';
    sd.insert_surf(surf);
  end
end

if ~doa_method_flag
  % Beamforming method
  try
    surf = sd.get_surf('bottom gt');
    if strcmpi(param.tomo_collate.surfData_mode,'overwrite')
      surf.y = NaN(size(twtt_bin));
      surf.y(ceil(Nsv/2)+1,:) = interp1(mdata.Time,1:length(mdata.Time),Bottom);
      sd.set_surf(surf);
    end
  catch ME
    surf = tomo.surfdata.empty_surf();
    surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
    surf.y = NaN(size(twtt_bin));
    surf.y(ceil(Nsv/2)+1,:) = interp1(mdata.Time,1:length(mdata.Time),Bottom);
    surf.plot_name_values = {'color','magenta','marker','+'};
    surf.name = 'bottom gt';
    sd.insert_surf(surf);
  end
else
  % DOA method
  try
    surf = sd.get_surf('bottom gt');
    if strcmpi(param.tomo_collate.surfData_mode,'overwrite')
      surf.y = NaN(size(twtt_bin));
      surf.y(1,:) = interp1(mdata.Time,1:length(mdata.Time),Bottom);
      sd.set_surf(surf);
    end
  catch ME
    surf = tomo.surfdata.empty_surf();
    surf.x = NaN(size(twtt_bin));
    for rline_idx = 1:Nx
      theta_rline = mdata.Tomo.theta(:,:,rline_idx);
      if ~all(isnan(theta_rline(:)))
        theta_rline = theta_rline(~isnan(theta_rline));
        surf.x(1:length(theta_rline),rline_idx) = theta_rline;
      end
    end
    surf.y = NaN(size(twtt_bin));
    surf.y(1,:) = interp1(mdata.Time,1:length(mdata.Time),Bottom);
    % Sort DOA min to max (and, accordingly, range-bins)
    [surf.x x_idx] = sort(surf.x*180/pi,1,'ascend');
    for rline_idx = 1:Nx
      surf.y(:,rline_idx) = surf.y(x_idx(:,rline_idx),rline_idx);
    end
  
    surf.plot_name_values = {'color','magenta','marker','+'};
    surf.name = 'bottom gt';
    sd.insert_surf(surf);
  end
end

if ~doa_method_flag
  % Beamforming method
  try
    surf = sd.get_surf('top gt');
    if strcmpi(param.tomo_collate.surfData_mode,'overwrite')
      surf.y = NaN(size(twtt_bin));
      surf.y(ceil(Nsv/2)+1,:) = interp1(mdata.Time,1:length(mdata.Time),Surface);
      sd.set_surf(surf);
    end
  catch ME
    surf = tomo.surfdata.empty_surf();
    surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
    surf.y = NaN * zeros(size(twtt_bin));
    surf.y(ceil(Nsv/2)+1,:) = interp1(mdata.Time,1:length(mdata.Time),Surface);
    surf.plot_name_values = {'color','magenta','marker','^'};
    surf.name = 'top gt';
    sd.insert_surf(surf);
  end
else
  try
    surf = sd.get_surf('top gt');
    if strcmpi(param.tomo_collate.surfData_mode,'overwrite')
      surf.y = NaN(size(twtt_bin));
      surf.y(1,:) = interp1(mdata.Time,1:length(mdata.Time),Surface);
      sd.set_surf(surf);
    end
  catch ME
    surf = tomo.surfdata.empty_surf();
    surf.x = NaN(size(twtt_bin));
    for rline_idx = 1:Nx
      theta_rline = mdata.Tomo.theta(:,:,rline_idx);
      if ~all(isnan(theta_rline(:)))
        theta_rline = theta_rline(~isnan(theta_rline));
        surf.x(1:length(theta_rline),rline_idx) = theta_rline;
      end
    end
    surf.y = NaN(size(twtt_bin));
    surf.y(1,:) = interp1(mdata.Time,1:length(mdata.Time),Surface);
    % Sort DOA min to max (and, accordingly, range-bins)
    [surf.x x_idx] = sort(surf.x*180/pi,1,'ascend');
    for rline_idx = 1:Nx
      surf.y(:,rline_idx) = surf.y(x_idx(:,rline_idx),rline_idx);
    end
  
    surf.plot_name_values = {'color','magenta','marker','^'};
    surf.name = 'top gt';
    sd.insert_surf(surf);
  end
end

if ~doa_method_flag
  % Beamforming method
  try
    surf = sd.get_surf('top quality');
    if strcmpi(param.tomo_collate.surfData_mode,'overwrite')
      surf.y = ones(size(twtt_bin));
      sd.set_surf(surf);
    end
  catch ME
    surf = tomo.surfdata.empty_surf();
    surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
    surf.y = ones(size(twtt_bin));
    surf.plot_name_values = {'color','red','marker','x'};
    surf.name = 'top quality';
    sd.insert_surf(surf);
  end
else
  try
    surf = sd.get_surf('top quality');
    if strcmpi(param.tomo_collate.surfData_mode,'overwrite')
      surf.y = NaN(size(twtt_bin));
%       surf.y = ones(size(twtt_bin));
      sd.set_surf(surf);
    end
  catch ME
    surf = tomo.surfdata.empty_surf();
    surf.x = NaN(size(twtt_bin));
    for rline_idx = 1:Nx
      theta_rline = mdata.Tomo.theta(:,:,rline_idx);
      if ~all(isnan(theta_rline(:)))
        theta_rline = theta_rline(~isnan(theta_rline));
        surf.x(1:length(theta_rline),rline_idx) = theta_rline;
      end
    end
%     surf.y = NaN(size(twtt_bin));
    surf.y = ones(size(twtt_bin));
    % Sort DOA min to max (and, accordingly, range-bins)
    [surf.x x_idx] = sort(surf.x*180/pi,1,'ascend');
    for rline_idx = 1:Nx
      surf.y(:,rline_idx) = surf.y(x_idx(:,rline_idx),rline_idx);
    end
  
    surf.plot_name_values = {'color','red','marker','x'};
    surf.name = 'top quality';
    sd.insert_surf(surf);
  end
end

if ~doa_method_flag
  % Beamforming method
  try
    surf = sd.get_surf('bottom quality');
    if strcmpi(param.tomo_collate.surfData_mode,'overwrite')
      surf.y = ones(size(twtt_bin));
      sd.set_surf(surf);
    end
  catch ME
    surf = tomo.surfdata.empty_surf();
    surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
    surf.y = ones(size(twtt_bin));
    surf.plot_name_values = {'color','red','marker','^'};
    surf.name = 'bottom quality';
    sd.insert_surf(surf);
  end
else
  try
    surf = sd.get_surf('bottom quality');
    if strcmpi(param.tomo_collate.surfData_mode,'overwrite')
      surf.y = ones(size(twtt_bin));
      sd.set_surf(surf);
    end
  catch ME
    surf = tomo.surfdata.empty_surf();
    surf.x = NaN(size(twtt_bin));
    for rline_idx = 1:Nx
      theta_rline = mdata.Tomo.theta(:,:,rline_idx);
      if ~all(isnan(theta_rline(:)))
        theta_rline = theta_rline(~isnan(theta_rline));
        surf.x(1:length(theta_rline),rline_idx) = theta_rline;
      end
    end
%     surf.y = NaN(size(twtt_bin));
     surf.y = ones(size(twtt_bin));
    % Sort DOA min to max (and, accordingly, range-bins)
    [surf.x x_idx] = sort(surf.x*180/pi,1,'ascend');
    for rline_idx = 1:Nx
      surf.y(:,rline_idx) = surf.y(x_idx(:,rline_idx),rline_idx);
    end
    
    surf.plot_name_values = {'color','red','marker','^'};
    surf.name = 'bottom quality';
    sd.insert_surf(surf);
  end
end

sd.set({'bottom','ice mask','bottom gt','bottom quality'}, ...
  'top','top','active','bottom','mask','ice mask','gt','bottom gt','quality','bottom quality');

sd.set({'top','top gt','top quality'}, ...
  'active','top','gt','top gt','quality','top quality');

sd.save_surfdata(out_fn,doa_method_flag);

mu = [];
sigma = [];

for cmd_idx = 1:length(param.tomo_collate.surfdata_cmds)
  cmd = param.tomo_collate.surfdata_cmds(cmd_idx).cmd;
  surf_names = param.tomo_collate.surfdata_cmds(cmd_idx).surf_names;
  if ischar(surf_names)
    surf_names = {surf_names};
  end
  visible = param.tomo_collate.surfdata_cmds(cmd_idx).visible;
  if isfield(param.tomo_collate.surfdata_cmds(cmd_idx),'plot_name_values') ...
      && ~isempty(param.tomo_collate.surfdata_cmds(cmd_idx).plot_name_values)
    plot_name_values = param.tomo_collate.surfdata_cmds(cmd_idx).plot_name_values;
  elseif exist('plot_name_values','var')
    % For DOA bottom tracking -- Do nothing
  else
    plot_name_values = {'color','black','marker','^'};
  end
  
  if any(strcmpi(cmd,{'detect','extract'})) && isempty(mu)
    %% Training parameters for image template's mu and sigma
    
    if isfield(param.tomo_collate.surfdata_cmds(cmd_idx),'data_threshold') ...
        && ~isempty(param.tomo_collate.surfdata_cmds(cmd_idx).data_threshold)
      data_threshold = param.tomo_collate.surfdata_cmds(cmd_idx).data_threshold;
    else
      data_threshold = 13.5;
    end
    
    for rline = 1:size(mdata.Tomo.img,3)
      if ~mod(rline-1,100)
        fprintf('  Training %d of %d (%s)\n', rline, size(mdata.Tomo.img,3), datestr(now));
      end
      thresh_data = data(:,:,rline);
      thresh_data(thresh_data>data_threshold) = data_threshold;
      [m, s] = tomo.train_model(thresh_data, ...
        double(twtt_bin(:,rline)));
      mu(:,rline) = m;
      sigma(:,rline) = s;
    end
    %   mu = repmat([11.6734184113247 11.8357634315107 11.7477015213467 11.5642270542054 11.3655718245298 11.2178010788707 11.11172116154 11.0442549382899 10.9800832668574 10.9047999009164 10.8000063888223],size(mdata.Tomo.img,3),1);
    %   sigma = repmat([17.9297614680615 18.5178215941504 17.1485050463076 15.8106528912151 14.7936777080171 14.146975962117 13.9673485950864 13.9574519525412 13.5837122364561 13.0310380580007 12.2855990897649],size(mdata.Tomo.img,3),1);
    mdata.Tomo.mu = mu;
    mdata.Tomo.sigma = sigma;
    
    % in_dir: Directory where 3D image files are at
    in_dir = ct_filename_out(param,param.tomo_collate.in_dir);
    
    % combined_fn: Filename with 3D data
    combined_fn = fullfile(in_dir,sprintf('Data_%s_%03.0f.mat',param.day_seg,param.load.frm));
    
    Tomo = mdata.Tomo;
    save(combined_fn,'-append','Tomo');
  end
  
  if strcmpi(cmd,'detect')
    %% Run detect
    
    if isfield(param.tomo_collate.surfdata_cmds(cmd_idx),'data_threshold') ...
        && ~isempty(param.tomo_collate.surfdata_cmds(cmd_idx).data_threshold)
      data_threshold = param.tomo_collate.surfdata_cmds(cmd_idx).data_threshold;
    else
      data_threshold = 13.5;
    end
    
    detect_surface = zeros(size(mdata.Tomo.img,2),size(mdata.Tomo.img,3));
    mu            = mdata.Tomo.mu;
    sigma         = mdata.Tomo.sigma;
    
    for rline = 1:size(mdata.Tomo.img,3)
      if ~mod(rline-1,500)
        fprintf('  Detect %d of %d (%s)\n', rline, size(mdata.Tomo.img,3), datestr(now));
      end
      
      % data_threshold: log scale data will be clipped to this threshold
      thresh_data = data(:,:,rline);
      thresh_data(thresh_data>data_threshold) = data_threshold;
      
      labels = tomo.detect(thresh_data, double(twtt_bin(:,rline)), double(Bottom_bin(rline)), ...
        [], double(ice_mask(:,rline)), double(mu(:,rline)), double(sigma(:,rline)));
      
      detect_surface(:,rline) = labels;
    end
    
    for surf_name_idx = 1:length(surf_names)
      surf_name = surf_names{surf_name_idx};
      try
        surf = sd.get_surf(surf_name);
        if ~strcmpi(param.tomo_collate.surfData_mode,'fillgaps')
          surf.y = detect_surface;
          surf.plot_name_values = plot_name_values;
          surf.visible = visible;
          sd.set_surf(surf);
        end
      catch ME
        surf = tomo.surfdata.empty_surf();
        surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
        surf.y = detect_surface;
        surf.name = surf_name;
        surf.plot_name_values = plot_name_values;
        surf.visible = visible;
        sd.insert_surf(surf);
      end
      sd.set(surf_name,'top','top','active','bottom','mask','ice mask', ...
        'gt','bottom gt','quality','bottom quality');
    end
    
  elseif strcmpi(cmd,'dem')
    %% Run DEM
    
    param.tomo_collate.surfdata_cmds(cmd_idx).dem_bad_value = -32767;
    param.tomo_collate.surfdata_cmds(cmd_idx).dem_guard = 12e3;
    param.tomo_collate.surfdata_cmds(cmd_idx).dem_per_slice_guard = 240;
    
    % Load Geotiff projection information
    proj = geotiffinfo(param.tomo_collate.surfdata_cmds(cmd_idx).dem_fn);
    
    % Load Geotiff raster, convert to double, set bad values to NaN
    [DEM, R, ~] = geotiffread(param.tomo_collate.surfdata_cmds(cmd_idx).dem_fn);
    DEM = double(DEM);
    DEM(DEM==param.tomo_collate.surfdata_cmds(cmd_idx).dem_bad_value) = NaN;
    
    % Create Geotiff axes
    DEM_x = R(3,1) + R(2,1)*(0:size(DEM,2)-1);
    DEM_y = R(3,2) + R(1,2)*(0:size(DEM,1)-1);
    
    DEM_x_mesh = repmat(DEM_x,[size(DEM,1) 1]);
    DEM_y_mesh= repmat(DEM_y',[1 size(DEM,2)]);
    
    % Project flight line to DEM coordinates
    [mdata.x,mdata.y] = projfwd(proj,mdata.Latitude,mdata.Longitude);
    
    % Process inputs
    dem_guard = param.tomo_collate.surfdata_cmds(cmd_idx).dem_guard;
    dem_per_slice_guard = param.tomo_collate.surfdata_cmds(cmd_idx).dem_per_slice_guard;
    
    top_idx = sd.get_index('top');
    
    theta = mdata.Tomo.theta(:,1);
    if isfield(param.tomo_collate,'sv_cal_fn') && ~isempty(param.tomo_collate.sv_cal_fn)
      theta_cal = load(param.tomo_collate.sv_cal_fn);
      theta = theta_cal.theta;
      theta_cal = theta;
    end
    
    % Loop to create dem_surface from DEM one range line at a time
    dem_surface = NaN*zeros(size(mdata.Tomo.img,2),size(mdata.Tomo.img,3));
    for rline = 1:size(mdata.Tomo.img,3)
      if ~mod(rline-1,500)
        fprintf('  DEM %d of %d (%s)\n', rline, size(mdata.Tomo.img,3), datestr(now));
      end
      
      DEM_mask = DEM_x_mesh > mdata.x(rline)-dem_guard & DEM_x_mesh < mdata.x(rline)+dem_guard ...
        & DEM_y_mesh > mdata.y(rline)-dem_guard & DEM_y_mesh < mdata.y(rline)+dem_guard ...
        & ~isnan(DEM);
      DEM_idxs = find(DEM_mask);
      
      if numel(DEM_idxs)==0
        % warning('Range Line %d of Frame %d is not spanned by DEM.',rline,param.load.frm);
        continue;
      end
      
      % Convert from projection to geodetic (lat,lon,elev)
      [DEM_lat,DEM_lon] = projinv(proj,DEM_x_mesh(DEM_idxs),DEM_y_mesh(DEM_idxs));
      DEM_elev = DEM(DEM_idxs);
      
      if all(isnan(DEM(DEM_idxs)))
        continue;
      end
      
      % Convert from geodetic (lat,lon,elev) to ECEF (x,y,z)
      physical_constants;
      [DEM_ecef_x,DEM_ecef_y,DEM_ecef_z] = geodetic2ecef(single(DEM_lat)/180*pi,single(DEM_lon)/180*pi,single(DEM_elev),WGS84.ellipsoid);
      
      origin = mdata.param_array.array_proc.fcs.origin(:,rline);
      
      % Convert from ECEF to FCS/SAR
      Tfcs_ecef = [mdata.param_array.array_proc.fcs.x(:,rline), ...
        mdata.param_array.array_proc.fcs.y(:,rline), ...
        mdata.param_array.array_proc.fcs.z(:,rline)];
      Tecef_fcs = inv(Tfcs_ecef);
      
      tmp = Tecef_fcs * [DEM_ecef_x.'-origin(1); DEM_ecef_y.'-origin(2); DEM_ecef_z.'-origin(3)];
      DEM_fcs_x = tmp(1,:);
      DEM_fcs_y = tmp(2,:);
      DEM_fcs_z = tmp(3,:);
      
      slice_mask = DEM_fcs_x > -dem_per_slice_guard & DEM_fcs_x < dem_per_slice_guard;
      
      x = DEM_fcs_x(slice_mask);
      y = DEM_fcs_y(slice_mask);
      z = DEM_fcs_z(slice_mask);
      
      if(numel(x)>=3)
        faces = delaunay(double(x),double(y));
        vertices = [double(x).' double(y).' double(z).'];  % vertices stored as Nx3 matrix
        vert1 = vertices(faces(:,1),:);
        vert2 = vertices(faces(:,2),:);
        vert3 = vertices(faces(:,3),:);
        
        twtt = NaN*zeros(size(mdata.Tomo.img,2),1);
        for theta_idx = 1:length(theta)
          % Find the point on the top surface
          top_twtt = interp1(1:length(mdata.Time),mdata.Time,sd.surf(top_idx).y(theta_idx,rline));
          top_orig = [0 sin(theta(theta_idx))*top_twtt*c/2 -cos(theta(theta_idx))*top_twtt*c/2];
          
          % Calculate refraction
          theta_refract = asin(sin(theta(theta_idx))/sqrt(er_ice));
          
          dir_v = [0 sin(theta_refract) -cos(theta_refract)];
          
          [Intersect, t] = TriangleRayIntersection(top_orig, dir_v, vert1, vert2, vert3);
          
          intersect_idx = find(Intersect);
          
          if isempty(intersect_idx)
            twtt(theta_idx) = NaN;
          else
            twtt(theta_idx) = top_twtt + t(intersect_idx(1))/(c/2/sqrt(er_ice));
          end
        end
        dem_surface(:,rline) = interp1(mdata.Time,1:length(mdata.Time),twtt);
        
      end
      
    end
    
    for surf_name_idx = 1:length(surf_names)
      surf_name = surf_names{surf_name_idx};
      try
        surf = sd.get_surf(surf_name);
        if ~strcmpi(param.tomo_collate.surfData_mode,'fillgaps')
          surf.y = dem_surface;
          surf.plot_name_values = plot_name_values;
          surf.visible = visible;
          sd.set_surf(surf);
        end
      catch ME
        surf = tomo.surfdata.empty_surf();
        surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
        surf.y = dem_surface;
        surf.name = surf_name;
        surf.plot_name_values = plot_name_values;
        surf.visible = visible;
        sd.insert_surf(surf);
      end
      sd.set(surf_name,'top','top','active','bottom','mask','ice mask', ...
        'gt','bottom gt','quality','bottom quality');
    end
    
  elseif strcmpi(cmd,'extract')
    %% Run extract
    fprintf('  Extract (%s)\n', datestr(now));
    
    if isfield(param.tomo_collate.surfdata_cmds(cmd_idx),'data_threshold') ...
        && ~isempty(param.tomo_collate.surfdata_cmds(cmd_idx).data_threshold)
      data_threshold = param.tomo_collate.surfdata_cmds(cmd_idx).data_threshold;
    else
      data_threshold = 13.5;
    end
    
    mu            = mdata.Tomo.mu;
    sigma         = mdata.Tomo.sigma;
    
    % data_threshold: log scale data will be clipped to this threshold
    thresh_data = data;
    thresh_data(thresh_data>data_threshold) = data_threshold;
    
    Extra_bin = [];
    
    extract_surface = tomo.extract(double(thresh_data), double(twtt_bin), double(Bottom_bin), ...
      double(Extra_bin), double(ice_mask), double(mean(mu,2)), double(mean(sigma,2)));
    
    extract_surface = reshape(extract_surface,size(mdata.Tomo.img,2),size(mdata.Tomo.img,3));   
    % Visualization of mean and variance vectors
    if 0
      figure; (plot(transition_mu)); hold on;
      plot(transition_sigma); xlim([1 64])
      legend('Mean', 'Variance', 'Location', 'northwest');
      xlabel('DoA bins');
    end
    
    
    for surf_name_idx = 1:length(surf_names)
      surf_name = surf_names{surf_name_idx};
      try
        surf = sd.get_surf(surf_name);
        if ~strcmpi(param.tomo_collate.surfData_mode,'fillgaps')
          surf.y = extract_surface;
          surf.plot_name_values = plot_name_values;
          surf.visible = visible;
          sd.set_surf(surf);
        end
      catch ME
        surf = tomo.surfdata.empty_surf();
        surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
        surf.y = extract_surface;
        surf.name = surf_name;
        surf.plot_name_values = plot_name_values;
        surf.visible = visible;
        sd.insert_surf(surf);
      end
      sd.set(surf_name,'top','top','active','bottom','mask','ice mask', ...
        'gt','bottom gt','quality','bottom quality');
    end
    
  elseif strcmpi(cmd,'viterbi')
    %% Run Viterbi
    viterbi_surface = zeros(size(mdata.Tomo.img,2),size(mdata.Tomo.img,3));
    % Check for smoothness weight
    if isfield(param.tomo_collate.surfdata_cmds(cmd_idx),'smooth_weight') ...
        && ~isempty(param.tomo_collate.surfdata_cmds(cmd_idx).smooth_weight)
      smooth_weight = param.tomo_collate.surfdata_cmds(cmd_idx).smooth_weight;
    else
      smooth_weight = 55; % schu
    end
    % Check for smoothness variance
    if isfield(param.tomo_collate.surfdata_cmds(cmd_idx),'smooth_var') ...
        && ~isempty(param.tomo_collate.surfdata_cmds(cmd_idx).smooth_var)
      smooth_var = param.tomo_collate.surfdata_cmds(cmd_idx).smooth_var;
    else
      smooth_var = -1;
    end
    % Check for repulsion weight
    if isfield(param.tomo_collate.surfdata_cmds(cmd_idx),'repulsion') ...
        && ~isempty(param.tomo_collate.surfdata_cmds(cmd_idx).repulsion)
      repulsion = param.tomo_collate.surfdata_cmds(cmd_idx).repulsion;
    else
      repulsion = 150; % schu
    end
    % Check for extra ground truth weight
    if isfield(param.tomo_collate.surfdata_cmds(cmd_idx),'egt_weight') ...
        && ~isempty(param.tomo_collate.surfdata_cmds(cmd_idx).egt_weight)
      egt_weight = param.tomo_collate.surfdata_cmds(cmd_idx).egt_weight;
    else
      egt_weight = 10;    
    % Visualization of mean and variance vectors
    if 0
      figure; (plot(transition_mu)); hold on;
      plot(transition_sigma); xlim([1 64])
      legend('Mean', 'Variance', 'Location', 'northwest');
      xlabel('DoA bins');
    end
    end
    % Check for ice_mask scanning threshold
    if isfield(param.tomo_collate.surfdata_cmds(cmd_idx),'ice_bin_thr') ...
        && ~isempty(param.tomo_collate.surfdata_cmds(cmd_idx).ice_bin_thr)
      ice_bin_thr = param.tomo_collate.surfdata_cmds(cmd_idx).ice_bin_thr;
    else
      ice_bin_thr = 3;
    end
    % Check for slope
    if isfield(param.tomo_collate.surfdata_cmds(cmd_idx),'slope') ...
        && ~isempty(param.tomo_collate.surfdata_cmds(cmd_idx).slope)
      slope = param.tomo_collate.surfdata_cmds(cmd_idx).slope;
    else
      slope = zeros(1, size(data,2)-1);
    end
    
    mu_size = 11;
    mu      = sinc(linspace(-1.5, 1.5, mu_size));
    sigma   = sum(mu)/20*ones(1,mu_size);
    
    %% Distance-to-Ice-Margin model
    clear DIM DIM_costmatrix;
    
    DIM = load(fullfile(param.path, '+tomo', 'Layer_tracking_3D_parameters_Matrix.mat'));
    DIM_costmatrix = DIM.Layer_tracking_3D_parameters;
    DIM_costmatrix = DIM_costmatrix .* (200 ./ max(DIM_costmatrix(:)));

    %% DoA-to-DoA transition model
    % Obtained from geostatistical analysis of 2014 Greenland P3
    transition_mu = [0.000000, 0.000000, 2.590611, 3.544282, 4.569263, 5.536577, 6.476430, 7.416807, 8.404554, 9.457255, 10.442658, 11.413710, 12.354409, 13.332689, 14.364614, 15.381671, 16.428969, 17.398906, 18.418794, 19.402757, 20.383026, 21.391834, 22.399259, 23.359765, 24.369957, 25.344982, 26.301805, 27.307530, 28.274756, 28.947572, 29.691010, 32.977387, 34.203212, 34.897994, 35.667128, 36.579019, 37.558978, 38.548659, 39.540715, 40.550138, 41.534781, 42.547407, 43.552700, 44.537758, 45.553618, 46.561057, 47.547331, 48.530976, 49.516588, 50.536075, 51.562886, 52.574938, 53.552979, 54.554206, 55.559657, 56.574029, 57.591999, 58.552986, 59.562937, 60.551616, 61.549909, 62.551092, 63.045791, 63.540490];
    transition_sigma = [0.457749, 0.805132, 1.152514, 1.213803, 1.290648, 1.370986, 1.586141, 1.626730, 1.785789, 1.791043, 1.782936, 1.727153, 1.770210, 1.714973, 1.687484, 1.663294, 1.633185, 1.647318, 1.619522, 1.626555, 1.649593, 1.628138, 1.699512, 1.749184, 1.809822, 1.946782, 2.126822, 2.237959, 2.313358, 2.280555, 1.419753, 1.112363, 1.426246, 2.159619, 2.140899, 2.083267, 1.687420, 1.574745, 1.480296, 1.443887, 1.415708, 1.356100, 1.401891, 1.398477, 1.365730, 1.418647, 1.407810, 1.430151, 1.391357, 1.403471, 1.454194, 1.470535, 1.417235, 1.455086, 1.436509, 1.378037, 1.415834, 1.333177, 1.298108, 1.277559, 1.358260, 1.483521, 1.674642, 1.865764];
    
    for rline = 1:size(mdata.Tomo.img,3)
      if ~mod(rline-1,500)
        fprintf('  Viterbi %d of %d (%s)\n', rline, size(mdata.Tomo.img,3), datestr(now));
      end
      
      detect_data    = data(:,:,rline);
      surf_bins      = twtt_bin(:,rline).';
      bottom_bin     = Bottom_bin(rline);
      gt             = [32; bottom_bin + 0.5];
      
      % Check for viterbi weight
      if isfield(param.tomo_collate.surfdata_cmds(cmd_idx),'viterbi_weight') ...
          && ~isempty(param.tomo_collate.surfdata_cmds(cmd_idx).viterbi_weight)
        viterbi_weight = param.tomo_collate.surfdata_cmds(cmd_idx).viterbi_weight;
      else
        viterbi_weight = ones([1 size(data,2)]);
        viterbi_weight(gt(1, :)) = 2;
      end
      
      threshold = 13.5;
      detect_data(detect_data>threshold) = threshold;
      detect_data = fir_dec(detect_data.',hanning(3).'/3,1).';
      
      mask = ice_mask(:,rline);
      mask_dist = round(bwdist(mask == 0));
      
      bounds = [1 length(surf_bins)];
      
      %% Call viterbi.cpp
      labels = tomo.viterbi(double(detect_data), double(surf_bins), ...
        double(bottom_bin), double(gt), double(mask), double(mu), ...
        double(sigma), double(egt_weight), double(smooth_weight), ...
        double(smooth_var), double(slope), int64(bounds), ...
        double(viterbi_weight), double(repulsion), double(ice_bin_thr), ...
        double(mask_dist), double(DIM_costmatrix), ...
        double(transition_mu), double(transition_sigma));
      
      viterbi_surface(:,rline) = labels;
    end
    
    for surf_name_idx = 1:length(surf_names)
      surf_name = surf_names{surf_name_idx};
      try
        surf = sd.get_surf(surf_name);
        if ~strcmpi(param.tomo_collate.surfData_mode,'fillgaps')
          surf.y = viterbi_surface;
          surf.plot_name_values = plot_name_values;
          surf.visible = visible;
          sd.set_surf(surf);
        end
      catch ME
        surf = tomo.surfdata.empty_surf();
        surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
        surf.y = viterbi_surface;
        surf.name = surf_name;
        surf.plot_name_values = plot_name_values;
        surf.visible = visible;
        sd.insert_surf(surf);
      end
      sd.set(surf_name,'top','top','active','bottom','mask','ice mask', ...
        'gt','bottom gt','quality','bottom quality');
    end
    
  elseif strcmpi(cmd,'trws')
    %% Run TRW-S
    fprintf('  TRW-S (%s)\n', datestr(now));
    % Check for smoothness weight
    if isfield(param.tomo_collate.surfdata_cmds(cmd_idx),'smooth_weight') ...
        && ~isempty(param.tomo_collate.surfdata_cmds(cmd_idx).smooth_weight)
      smooth_weight = param.tomo_collate.surfdata_cmds(cmd_idx).smooth_weight;
    else
      smooth_weight = [22 22];
    end
    % Check for smoothness variance
    if isfield(param.tomo_collate.surfdata_cmds(cmd_idx),'smooth_var') ...
        && ~isempty(param.tomo_collate.surfdata_cmds(cmd_idx).smooth_var)
      smooth_var = param.tomo_collate.surfdata_cmds(cmd_idx).smooth_var;
    else
      smooth_var = 32;
    end
    % Check for max number of loops
    if isfield(param.tomo_collate.surfdata_cmds(cmd_idx),'max_loops') ...
        && ~isempty(param.tomo_collate.surfdata_cmds(cmd_idx).max_loops)
      max_loops = param.tomo_collate.surfdata_cmds(cmd_idx).max_loops;
    else
      max_loops = 50;
    end
    
    smooth_slope = [];
    mu_size = 11;
    mu = sinc(linspace(-1.5,1.5,mu_size));
    sigma = sum(mu)/20*ones(1,mu_size);
    bounds = [param.tomo_collate.bounds_relative(1) size(data,2)-1-param.tomo_collate.bounds_relative(2) -1 -1];    
    %mask_dist = round(bwdist(ice_mask == 0));
    mask_dist = inf*ones(size(ice_mask));
    mask_dist = round(mask_dist .* 9);
    clear DIM DIM_costmatrix;
    
    DIM = load(fullfile(param.path, '+tomo', 'Layer_tracking_3D_parameters_Matrix.mat'));
    DIM_costmatrix = DIM.Layer_tracking_3D_parameters;
    DIM_costmatrix = DIM_costmatrix .* (200 ./ max(DIM_costmatrix(:)));

    
    %% DoA-to-DoA transition model
    % Obtained from geostatistical analysis of 2014 Greenland P3
    transition_mu = [2.0436 2.3331 2.5009 3.3719 4.6784 5.6978 6.5621 7.5174 8.5156 9.5651 10.5363 11.5323 12.5066 13.5002 14.4998 15.5585 16.5564 17.5435 18.5288 19.5175 20.5071 21.5108 22.5106 23.4993   24.4847 25.4574 26.4393 27.4864 28.4248 29.1076 29.7335 32.9690 34.1460 34.6690 35.4782 36.4208 37.4689 38.4754 39.4688 40.4474 41.4559 42.4452 43.4168 44.4374 45.4158 46.4087 47.4159 48.4306 49.4311 50.4148 51.4397 52.4642 53.4303 54.4758 55.4716 56.4896 57.5388 58.5285 59.4507 60.4436 61.4986 62.5633 62.6210 62.6788];
    transition_sigma = [1.9131 1.3669 1.5377 1.7085 1.7066 1.8079 1.8992 2.0351 2.0593 2.0194 1.9234 1.9222 1.9576 1.8838 1.9062 1.8439 1.7892 1.7756 1.7726 1.8337 1.7814 1.8196 1.9341 1.9805 2.1382 2.2869 2.4564 2.4599 2.4413 2.3801 1.4076 1.0751 1.3504 1.8570 2.0304 2.1111 1.8376 1.6472 1.5613 1.5116 1.4367 1.4435 1.4491 1.4410 1.4299 1.4022 1.4598 1.4219 1.4193 1.4158 1.4456 1.4779 1.4647 1.5021 1.4541 1.4040 1.4053 1.2808 1.2195 1.1342 1.3246 1.2063 1.6347 2.0632];
    RLINE_transition_sigma = [733.371814 77.126528 39.263353 12.295813 4.374837 6.958925 5.930228 2.258107 1.428613 1.388027 0.752566 0.979279 0.619339 0.763956 0.617092 0.627093 0.535119 0.488883 0.466207 0.452399 0.449242 0.448823 0.436793 0.434663 0.423571 0.446478 0.438191 0.429951 0.423105 0.403451 0.391212 0.375924 0.386908 0.385439 0.395184 0.401038 0.402252 0.409454 0.411048 0.414927 0.415572 0.421242 0.438985 0.455439 0.473058 0.490692 0.512094 0.544576 0.593388 0.629358 0.624657 0.648194 0.696526 0.754032 0.842585 0.960767 1.158873 1.759647 2.576886 5.313411 8.513457 18.789263 45.256746 82.139656];
    
    if length(transition_mu) ~= Nsv
      transition_mu = imresize(transition_mu, [1 Nsv]);
    end
    
    if length(transition_sigma) ~= Nsv
      transition_sigma = imresize(transition_sigma, [1 Nsv]);
    end
    
    if length(RLINE_transition_sigma) ~= Nsv
      RLINE_transition_sigma = imresize(RLINE_transition_sigma, [1 Nsv]);
    end
    % Visualization of mean and variance vectors
    if 0
      figure; (plot(transition_mu)); hold on;
      plot(transition_sigma); xlim([1 64])
      legend('Mean', 'Variance', 'Location', 'northwest');
      xlabel('DoA bins');
    end
    smooth_weight = 0.08 .* smooth_weight;
    gt = zeros(3, length(Bottom_bin));
    gt(1, :) = 1 : length(Bottom_bin);
    gt(2, :) = round(Nsv ./ 2) * ones(1, length(Bottom_bin));
    gt(3, :) = Bottom_bin(:) + 0.5;

    tic;
    trws_surface = tomo.trws(double(data), ...
      double(twtt_bin), double(Bottom_bin), double(gt), double(ice_mask), ...
      double(mu), double(sigma), double(smooth_weight), double(smooth_var), ...
      double(smooth_slope), double([]), double(max_loops), int64(bounds), ...
      double(mask_dist), double(DIM_costmatrix), ...
      double(transition_mu), double(transition_sigma), ...
      double(RLINE_transition_sigma));
    toc;
    
    trws_surface = reshape(trws_surface,size(mdata.Tomo.img,2), ...
      size(mdata.Tomo.img,3));
    
    for surf_name_idx = 1:length(surf_names)
      surf_name = surf_names{surf_name_idx};
      try
        surf = sd.get_surf(surf_name);
        if ~strcmpi(param.tomo_collate.surfData_mode,'fillgaps')
          surf.y = trws_surface;
          surf.plot_name_values = plot_name_values;
          surf.visible = visible;
          sd.set_surf(surf);
        end
      catch ME
        surf = tomo.surfdata.empty_surf();
        surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
        surf.y = trws_surface;
        surf.name = surf_name;
        surf.plot_name_values = plot_name_values;
        surf.visible = visible;
        sd.insert_surf(surf);
      end
      sd.set(surf_name,'top','top','active','bottom','mask','ice mask', ...
        'gt','bottom gt','quality','bottom quality');
    end
    
  elseif strcmpi(cmd,'c3d_rnn')
    c3d_rnn.dwnsammat_dir = fullfile(ct_filename_out(param, 'C3D_RNN_temporary_resources'), '');
    c3d_rnn.dwnsamnpy_dir = fullfile(ct_filename_out(param, 'C3D_RNN_temporary_resources'), '');
    temp_str              = strfind(c3d_rnn.dwnsammat_dir, filesep);
    c3d_rnn.dwnsammat_dir = c3d_rnn.dwnsammat_dir(1 : temp_str(end));
    c3d_rnn.dwnsammat_dir = fullfile(c3d_rnn.dwnsammat_dir, 'slices_mat_64x64', ct_filename_out(param, filesep));
    c3d_rnn.dwnsamnpy_dir = c3d_rnn.dwnsamnpy_dir(1 : temp_str(end));
    c3d_rnn.dwnsamnpy_dir = fullfile(c3d_rnn.dwnsamnpy_dir, 'slices_npy_64x64', ct_filename_out(param, filesep));
    
    % Down-sample data and save file for each and every slice
    for rline = 1 : size(data, 3)
      fusion = db(mdata.Topography.img(:, :, rline));
      fusion(fusion>27) = 27;
      fusion = imresize(fusion, [64, 64]);
      fusion = mat2gray(fusion);
      
      out_dir = char(strcat(c3d_rnn.dwnsammat_dir, '/', sprintf('%.3d', param.proc.frm), '/'));
      outfile = strcat(out_dir, num2str(rline,'%05.f'), '.mat');
      
      try
        save(outfile, 'fusion');
      catch ME
        try
          mkdir(out_dir);
          save(outfile, 'fusion');
        catch ME
          fprintf('\nProblem saving temporary MAT file, verify.\n');
          keyboard
        end
      end
    end
    fprintf('\nFinished down-sampling and saving MAT files for %s_%03.0f.\nPath: %s\n\n',param.day_seg,param.proc.frm, out_dir);
    
    %% Search for pre-trained model files (c3d.pth and rnn.pth)
    c3d_rnn.pth_path     = fullfile(param.path, '+tomo', 'c3d_rnn_models');
    c3d_rnn.c3d_pth_path = fullfile(c3d_rnn.pth_path, 'c3d.pth');
    c3d_rnn.rnn_pth_path = fullfile(c3d_rnn.pth_path, 'rnn.pth');
    
    if ~exist(c3d_rnn.c3d_pth_path, 'file') || ~exist(c3d_rnn.c3d_pth_path, 'file')
      fprintf('\nProblem locating model files (c3d.pth and rnn.pth). Directories being searched:\n%s\n%s\n', ...
        c3d_rnn.c3d_pth_path, c3d_rnn.rnn_pth_path);
      keyboard
    end
    
    %% Convert from MAT to NPY and run C3D_RNN
    %   Calls shell script
    fprintf('Executing shell script to run Python scripts...\n\n');
    try
      temp_str              = strfind(c3d_rnn.dwnsammat_dir, filesep);
      c3d_rnn.dwnsammat_dir = c3d_rnn.dwnsammat_dir(1 : temp_str(end));
      c3d_rnn.dwnsamnpy_dir = c3d_rnn.dwnsamnpy_dir(1 : temp_str(end));
      temp_str              = strfind(c3d_rnn.dwnsammat_dir, filesep);
      c3d_rnn.data_dir      = c3d_rnn.dwnsammat_dir(1 : temp_str(end - 1));
      c3d_rnn.outtext_dir   = fullfile(c3d_rnn.data_dir, sprintf('C3D_RNN_Out_%s_%03.0f.txt', ...
        param.day_seg,param.proc.frm));
      cd(fullfile(param.path, '+tomo'));
      sh_cmd = ['chmod +x run_C3D_RNN_Python.sh; ./run_C3D_RNN_Python.sh ', ...
        c3d_rnn.dwnsammat_dir, ' ', ...
        c3d_rnn.dwnsamnpy_dir, ' ', ...
        c3d_rnn.data_dir,      ' ', ...
        c3d_rnn.c3d_pth_path,  ' ', ...
        c3d_rnn.rnn_pth_path,  ' ', ...
        c3d_rnn.outtext_dir,   ' ', ...
        sprintf('Data_%s_%03.0f',param.day_seg,param.proc.frm)];
      tic;
      system(sh_cmd);
      toc;
    catch ME
      fprintf('\nProblem during execution of Python scripts, verify.\n');
      keyboard
    end
    fprintf('\nFinished executing Python scripts for %s_%03.0f.\n\n',param.day_seg,param.proc.frm);
    
    c3d_rnn.result_surface = ones(size(mdata.Topography.img,2), size(mdata.Topography.img,3));
    c3d_rnn.result_bottom  = ones(size(mdata.Topography.img,2), size(mdata.Topography.img,3));
    
    sl_idx = 1;
    %% Get surface and bottom vectors from generated text file
    fid = fopen(c3d_rnn.outtext_dir);
    tline = '';
    
    while ischar(tline)
      tline = fgetl(fid);
      if (isscalar(tline) && tline == -1) || isempty(tline)
        continue;
      end
      tline = fgetl(fid);
      if (isscalar(tline) && tline == -1) || isempty(tline)
        continue;
      end
      c3d_rnn.result_surface(:, sl_idx) = str2num(tline)';
      tline = fgetl(fid);
      if (isscalar(tline) && tline == -1) || isempty(tline)
        continue;
      end
      c3d_rnn.result_bottom(:, sl_idx) = str2num(tline)';
      sl_idx = sl_idx + 1;
    end
    
    % Threshold bottom layer to location of surface layer (avoids negative
    % ice thickness)
    if param.tomo_collate.surfdata_cmds(cmd_idx).surface_threshold
      c3d_rnn.result_bottom(c3d_rnn.result_bottom < c3d_rnn.result_surface) = ...
        c3d_rnn.result_surface(c3d_rnn.result_bottom < c3d_rnn.result_surface);
    end
    
    c3d_rnn.result_matrix(:, :, 1) = c3d_rnn.result_surface;
    c3d_rnn.result_matrix(:, :, 2) = c3d_rnn.result_bottom;
    
    for surf_name_idx = 1:length(surf_names)
      surf_name = surf_names{surf_name_idx};
      try
        surf = sd.get_surf(surf_name);
        if ~strcmpi(param.tomo_collate.surfData_mode,'fillgaps')
          surf.y = c3d_rnn.result_matrix(:, :, surf_name_idx);
          surf.plot_name_values = plot_name_values;
          surf.visible = visible;
          sd.set_surf(surf);
        end
      catch ME
        surf = tomo.surfdata.empty_surf();
        surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
        surf.y = c3d_rnn.result_matrix(:, :, surf_name_idx);
        surf.name = surf_name;
        surf.plot_name_values = plot_name_values;
        surf.visible = visible;
        sd.insert_surf(surf);
      end
      sd.set(surf_name,'top','top','active','bottom','mask','ice mask', ...
        'gt','bottom gt','quality','bottom quality');
    end
  elseif strcmpi(cmd,'doa')
    %% DOA method: Only 'bottom' is supported at this point
    doa_surface.x = NaN(size(twtt_bin));
    doa_surface.y = NaN(size(twtt_bin));
    for rline_idx = 1:Nx
      theta_rline = mdata.Tomo.theta(:,:,rline_idx);
      if ~all(isnan(theta_rline(:)))
        [theta_rline_r theta_rline_c] = find(~isnan(theta_rline));
        theta_rline = theta_rline(~isnan(theta_rline));
        % Sort DOA min to max (and, accordingly, range-bins)
%         doa_surface.x(1:length(theta_rline),rline_idx) = theta_rline;
        [doa_surface.x(1:length(theta_rline),rline_idx), x_idx] = sort(theta_rline*180/pi,'ascend');
        doa_surface.y(x_idx,rline_idx) = theta_rline_r;
      end
    end
    
    % Sort DOA min to max (and, accordingly, range-bins)
%     [doa_surface.x x_idx] = sort(doa_surface.x*180/pi,1,'ascend');
%     for rline_idx = 1:Nx
%       doa_surface.y(:,rline_idx) = doa_surface.y(x_idx(:,rline_idx),rline_idx);
%     end
    
    for surf_name_idx = 1:length(surf_names)
      surf_name = surf_names{surf_name_idx};
      try
        surf = sd.get_surf(surf_name);
        if ~strcmpi(param.tomo_collate.surfData_mode,'fillgaps')
          surf.y = doa_surface.y;
          % Ensure non-negative ice thickness
          if merge_bottom_above_top && exist('ice_top','var') && isfield(ice_top,'y') && ~isempty(ice_top.y)
            surf.y(surf.y<ice_top.y) = ice_top.y(surf.y<ice_top.y);
          end
          surf.plot_name_values = plot_name_values;
          surf.visible = visible;
          sd.set_surf(surf);
        end
      catch ME
        surf = tomo.surfdata.empty_surf();
        surf.x = doa_surface.x;
        surf.y = doa_surface.y;
        % Ensure non-negative ice thickness
        if merge_bottom_above_top && exist('ice_top','var') && isfield(ice_top,'y') && ~isempty(ice_top.y)
          surf.y(surf.y<ice_top.y) = ice_top.y(surf.y<ice_top.y);
        end
        surf.name = surf_name;
        surf.plot_name_values = plot_name_values;
        surf.visible = visible;
        sd.insert_surf(surf);
      end
      sd.set(surf_name,'top','top','active','bottom','mask','ice mask', ...
        'gt','bottom gt','quality','bottom quality');
    end
  end
end
% end

sd.save_surfdata(out_fn,doa_method_flag);
fprintf('Done (%s)\n', datestr(now));


