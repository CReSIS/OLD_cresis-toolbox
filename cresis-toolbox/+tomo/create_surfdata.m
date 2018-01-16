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
%     .cmd: String with command: 'detect', 'extract', 'viterbi', or 'trws'
%     .surf_names: String or cell array of strings with surface names to be
%       updated by the results of the command
%     .visible: The visibility setting for this layer
%     DETECT parameters
%     .data_threshold: pixel values above this will be clipped (default is 13.5)
%     EXTRACT parameters
%     .data_threshold: pixel values above this will be clipped (default is 13.5)
%     TRWS parameters
%     .smooth_weight: smoothness weight [] (default is [22 22])
%     .smooth_var: Gaussian weighting in the elevation angle bin dimension (default is 32)
%     .max_loops: number of loops to run for (default is 50)
%   .bounds_relative: DOA bins and along-track slices to trim off from each edge [top bottom left right]
%     E.g. [0 0 0 0] would use all pixels, [1 0 0 0] would ignore the first row.
%     Usually trimming a few rows off the top/bottom is a good idea if
%     these represent angles out to +/- 90 deg.
%
% mdata: 3D data file struct (Data_YYYYMMDD_SS_FFF.mat with Topography
%   field)
%
% Outputs:
%   NONE
%
% See also: tomo.run_collate, tomo.collate, tomo_collate_task,
%   tomo.fuse_images, tomo.add_icemask_surfacedem, tomo.create_surfdata,
%
% Author: John Paden, Jordan Sprick, and Mingze Xu

if ~isfield(param.tomo_collate,'out_dir') || isempty(param.tomo_collate.out_dir)
  param.tomo_collate.out_dir = 'surfData';
end

%% Load surface and bottom information
param_load_layers = param;
param_load_layers.cmd.frms = round([-1,0,1] + param.proc.frm);

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

if length(Bottom)~=size(mdata.Topography.img,3)
  error('This should not happen.');
  Bottom = mdata.Bottom;
end
if length(Surface)~=size(mdata.Topography.img,3)
  error('This should not happen.');
  Surface = mdata.Surface;
end

%% Interpolate Bottom, mdata.twtt from twtt to bins
Bottom_bin = interp1(mdata.Time, 1:length(mdata.Time), Bottom);
if isfield(mdata,'ice_mask')
  ice_mask = mdata.ice_mask;
else
  ice_mask = ones(size(mdata.twtt));
end
Bottom_bin(isnan(Bottom_bin)) = -1;

%% Surface tracking prep: Convert img to double and log-scale
data = 10*log10(double(mdata.Topography.img));

%% Surface tracking prep
% 1. Convert from twtt to bins
twtt_bin = round(interp1(mdata.Time, 1:length(mdata.Time), mdata.twtt));
% 2. The tracking software assumes that the surface never approaches
%    within mu_lengthof the top/bottom of the range line, so we truncate
%    surface to ensure this never happens.
mu_length = 11;
twtt_bin(isnan(twtt_bin) | twtt_bin > length(mdata.Time)-mu_length) = length(mdata.Time)-mu_length;

%% Create output filename
out_dir = ct_filename_out(param,param.tomo_collate.out_dir,'');
if ~isdir(out_dir)
  mkdir(out_dir);
end
out_fn_name = sprintf('Data_%s_%03.0f.mat',param.day_seg,param.proc.frm);
out_fn = fullfile(out_dir,out_fn_name);
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
  sd.radar_name = mdata.param_combine.radar_name;
  sd.season_name = mdata.param_combine.season_name;
  sd.day_seg = mdata.param_combine.day_seg;
  sd.frm = mdata.param_combine.combine.frm;
  sd.gps_time = mdata.GPS_time;
  sd.theta = mdata.theta(:); % Make a column vector
  sd.time = mdata.Time(:); % Make a column vector
  sd.FCS.origin = mdata.param_combine.array_param.fcs{1}{1}.origin;
  sd.FCS.x = mdata.param_combine.array_param.fcs{1}{1}.x;
  sd.FCS.y = mdata.param_combine.array_param.fcs{1}{1}.y;
  sd.FCS.z = mdata.param_combine.array_param.fcs{1}{1}.z;
end

Nsv = size(mdata.Topography.img,2);

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

try
  surf = sd.get_surf('bottom');
  if strcmpi(param.tomo_collate.surfData_mode,'overwrite')
    surf.y = NaN * zeros(size(twtt_bin));
    sd.set_surf(surf);
  end
catch ME
  surf = tomo.surfdata.empty_surf();
  surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
  surf.y = NaN * zeros(size(twtt_bin));
  surf.plot_name_values = {'color','blue','marker','^'};
  surf.name = 'bottom';
  sd.insert_surf(surf);
end

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

try
  surf = sd.get_surf('bottom gt');
  if strcmpi(param.tomo_collate.surfData_mode,'overwrite')
    surf.y = NaN * zeros(size(twtt_bin));
    surf.y(ceil(Nsv/2)+1,:) = interp1(mdata.Time,1:length(mdata.Time),Bottom);
    sd.set_surf(surf);
  end
catch ME
  surf = tomo.surfdata.empty_surf();
  surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
  surf.y = NaN * zeros(size(twtt_bin));
  surf.y(ceil(Nsv/2)+1,:) = interp1(mdata.Time,1:length(mdata.Time),Bottom);
  surf.plot_name_values = {'color','magenta','marker','+'};
  surf.name = 'bottom gt';
  sd.insert_surf(surf);
end

try
  surf = sd.get_surf('top gt');
  if strcmpi(param.tomo_collate.surfData_mode,'overwrite')
    surf.y = NaN * zeros(size(twtt_bin));
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

sd.set({'bottom','ice mask','bottom gt','bottom quality'}, ...
  'top','top','active','bottom','mask','ice mask','gt','bottom gt','quality','bottom quality');

sd.set({'top','top gt','top quality'}, ...
  'active','top','gt','top gt','quality','top quality');

sd.save_surfdata(out_fn);

mu = [];
sigma = [];

for cmd_idx = 1:length(param.tomo_collate.surfdata_cmds)
  cmd = param.tomo_collate.surfdata_cmds(cmd_idx).cmd;
  surf_names = param.tomo_collate.surfdata_cmds(cmd_idx).surf_names;
  if ischar(surf_names)
    surf_names = {surf_names};
  end
  visible = param.tomo_collate.surfdata_cmds(cmd_idx).visible;
  
  if any(strcmpi(cmd,{'detect','extract'})) && isempty(mu)
    %% Training parameters for image template's mu and sigma
    
    if isfield(param.tomo_collate.surfdata_cmds(cmd_idx),'data_threshold') ...
        && ~isempty(param.tomo_collate.surfdata_cmds(cmd_idx).data_threshold)
      data_threshold = param.tomo_collate.surfdata_cmds(cmd_idx).data_threshold;
    else
      data_threshold = 13.5;
    end
    
    for rline = 1:size(mdata.Topography.img,3)
      if ~mod(rline-1,100)
        fprintf('  Training %d of %d (%s)\n', rline, size(mdata.Topography.img,3), datestr(now));
      end
      thresh_data = data(:,:,rline);
      thresh_data(thresh_data>data_threshold) = data_threshold;
      [m, s] = tomo.train_model(thresh_data, ...
        double(twtt_bin(:,rline)));
      mu(:,rline) = m;
      sigma(:,rline) = s;
    end
    %   mu = repmat([11.6734184113247 11.8357634315107 11.7477015213467 11.5642270542054 11.3655718245298 11.2178010788707 11.11172116154 11.0442549382899 10.9800832668574 10.9047999009164 10.8000063888223],size(mdata.Topography.img,3),1);
    %   sigma = repmat([17.9297614680615 18.5178215941504 17.1485050463076 15.8106528912151 14.7936777080171 14.146975962117 13.9673485950864 13.9574519525412 13.5837122364561 13.0310380580007 12.2855990897649],size(mdata.Topography.img,3),1);
    mdata.Topography.mu = mu;
    mdata.Topography.sigma = sigma;
    
    % in_dir: Directory where 3D image files are at
    in_dir = ct_filename_out(param,param.tomo_collate.in_dir);
    
    % combined_fn: Filename with 3D data
    combined_fn = fullfile(in_dir,sprintf('Data_%s_%03.0f.mat',param.day_seg,param.proc.frm));
    
    Topography = mdata.Topography;
    save(combined_fn,'-append','Topography');
  end
  
  if strcmpi(cmd,'detect')
    %% Run detect
    
    if isfield(param.tomo_collate.surfdata_cmds(cmd_idx),'data_threshold') ...
        && ~isempty(param.tomo_collate.surfdata_cmds(cmd_idx).data_threshold)
      data_threshold = param.tomo_collate.surfdata_cmds(cmd_idx).data_threshold;
    else
      data_threshold = 13.5;
    end
    
    detect_surface = zeros(size(mdata.Topography.img,2),size(mdata.Topography.img,3));
    mu            = mdata.Topography.mu;
    sigma         = mdata.Topography.sigma;
    
    for rline = 1:size(mdata.Topography.img,3)
      if ~mod(rline-1,500)
        fprintf('  Detect %d of %d (%s)\n', rline, size(mdata.Topography.img,3), datestr(now));
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
          sd.set_surf(surf);
        end
      catch ME
        surf = tomo.surfdata.empty_surf();
        surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
        surf.y = detect_surface;
        surf.plot_name_values = {'color','green','marker','^'};
        surf.name = surf_name;
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
    
    mu            = mdata.Topography.mu;
    sigma         = mdata.Topography.sigma;
    
    % data_threshold: log scale data will be clipped to this threshold
    thresh_data = data;
    thresh_data(thresh_data>data_threshold) = data_threshold;
    
    Extra_bin = [];
    
    extract_surface = tomo.extract(double(thresh_data), double(twtt_bin), double(Bottom_bin), ...
      double(Extra_bin), double(ice_mask), double(mean(mu,2)), double(mean(sigma,2)));
    
    extract_surface = reshape(extract_surface,size(mdata.Topography.img,2),size(mdata.Topography.img,3));
    
    for surf_name_idx = 1:length(surf_names)
      surf_name = surf_names{surf_name_idx};
      try
        surf = sd.get_surf(surf_name);
        if ~strcmpi(param.tomo_collate.surfData_mode,'fillgaps')
          surf.y = extract_surface;
          sd.set_surf(surf);
        end
      catch ME
        surf = tomo.surfdata.empty_surf();
        surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
        surf.y = extract_surface;
        surf.plot_name_values = {'color','yellow','marker','^'};
        surf.name = surf_name;
        surf.visible = visible;
        sd.insert_surf(surf);
      end
      sd.set(surf_name,'top','top','active','bottom','mask','ice mask', ...
        'gt','bottom gt','quality','bottom quality');
    end
    
  elseif strcmpi(cmd,'viterbi')
    %% Run Viterbi
    viterbi_surface = zeros(size(mdata.Topography.img,2),size(mdata.Topography.img,3));
    bounds = [param.tomo_collate.bounds_relative(1) size(data,2)-1-param.tomo_collate.bounds_relative(2)];
    for rline = 1:size(mdata.Topography.img,3)
      if ~mod(rline-1,500)
        fprintf('  Viterbi %d of %d (%s)\n', rline, size(mdata.Topography.img,3), datestr(now));
      end
      
      mu_size       = 11;
      mu            = sinc(linspace(-1.5, 1.5, mu_size));
      sigma         = sum(mu)/20*ones(1,mu_size);
      smooth_var    = -1;
      smooth_weight = 45; % 55
      repulsion     = 250; % 150
      smooth_weight = 55; % schu
      repulsion     = 150; % schu
      ice_bin_thr   = 3;
      
      detect_data = data(:,:,rline);
      surf_bins = twtt_bin(:,rline).';
      bottom_bin = Bottom_bin(rline);
      gt = [33; bottom_bin];
      viterbi_weight = ones([1 size(data,2)]);
      viterbi_weight(gt(1,:)) = 2;
      mask = ice_mask(:,rline).';
      egt_weight = 10;
      slope = zeros(1,63);
      
      labels = tomo.viterbi(double(detect_data), ...
        double(surf_bins), double(bottom_bin), ...
        double(gt), double(mask), ...
        double(mu), double(sigma), double(egt_weight), ...
        double(smooth_weight), double(smooth_var), double(slope), ...
        int64(bounds), double(viterbi_weight), ...
        double(repulsion), double(ice_bin_thr), 1);
      
      viterbi_surface(:,rline) = labels;
    end
    
    for surf_name_idx = 1:length(surf_names)
      surf_name = surf_names{surf_name_idx};
      try
        surf = sd.get_surf(surf_name);
        if ~strcmpi(param.tomo_collate.surfData_mode,'fillgaps')
          surf.y = viterbi_surface;
          sd.set_surf(surf);
        end
      catch ME
        surf = tomo.surfdata.empty_surf();
        surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
        surf.y = viterbi_surface;
        surf.plot_name_values = {'color','yellow','marker','^'};
        surf.name = surf_name;
        surf.visible = visible;
        sd.insert_surf(surf);
      end
      sd.set(surf_name,'top','top','active','bottom','mask','ice mask', ...
        'gt','bottom gt','quality','bottom quality');
    end
    
  elseif strcmpi(cmd,'trws')
    %% Run TRW-S
    fprintf('  TRW-S (%s)\n', datestr(now));
    
    if isfield(param.tomo_collate.surfdata_cmds(cmd_idx),'smooth_weight') ...
        && ~isempty(param.tomo_collate.surfdata_cmds(cmd_idx).smooth_weight)
      smooth_weight = param.tomo_collate.surfdata_cmds(cmd_idx).smooth_weight;
    else
      smooth_weight = [22 22];
    end
    if isfield(param.tomo_collate.surfdata_cmds(cmd_idx),'smooth_var') ...
        && ~isempty(param.tomo_collate.surfdata_cmds(cmd_idx).smooth_var)
      smooth_var = param.tomo_collate.surfdata_cmds(cmd_idx).smooth_var;
    else
      smooth_var = 32;
    end
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
    ice_mask_transition = 90*fir_dec(fir_dec(double(shrink(ice_mask,2)),ones(1,5)/3.7).',ones(1,5)/3.7).';
    ice_mask_transition(ice_mask_transition>=90) = inf;
    bounds = [param.tomo_collate.bounds_relative(1) size(data,2)-1-param.tomo_collate.bounds_relative(2) -1 -1];
    
    trws_surface = tomo.trws(data, double(twtt_bin), double(Bottom_bin), ...
      double([]), double(ice_mask_transition), double(mu), double(sigma), ...
      smooth_weight, smooth_var, double(smooth_slope), [], ...
      double(max_loops), int64(bounds));
    
    trws_surface = reshape(trws_surface,size(mdata.Topography.img,2),size(mdata.Topography.img,3));
    
    for surf_name_idx = 1:length(surf_names)
      surf_name = surf_names{surf_name_idx};
      try
        surf = sd.get_surf(surf_name);
        if ~strcmpi(param.tomo_collate.surfData_mode,'fillgaps')
          surf.y = trws_surface;
          sd.set_surf(surf);
        end
      catch ME
        surf = tomo.surfdata.empty_surf();
        surf.x = repmat((1:Nsv).',[1 size(mdata.twtt,2)]);
        surf.y = trws_surface;
        surf.plot_name_values = {'color','green','marker','^'};
        surf.name = surf_name;
        surf.visible = visible;
        sd.insert_surf(surf);
      end
      sd.set(surf_name,'top','top','active','bottom','mask','ice mask', ...
        'gt','bottom gt','quality','bottom quality');
    end
  end
end

fprintf('Done (%s)\n', datestr(now));

sd.save_surfdata(out_fn);
