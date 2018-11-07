
%% User Setup
% =====================================================================
param_override = [];

params = read_param_xls(ct_filename_param('snow_param_2011_Greenland_P3.xls'));

% Example to run a specific segment and frame by overriding parameter spreadsheet values
params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20110329_01');
% params = ct_set_params(params,'cmd.frms',[1]);

params = ct_set_params(params,'get_echogram_stats.data_type','CSARP_post/qlook');
params = ct_set_params(params,'get_echogram_stats.echogram_img',0);
params = ct_set_params(params,'get_echogram_stats.noise_bins',[-400 -100]);
params = ct_set_params(params,'get_echogram_stats.signal_bins',[-99 500]);

dbstop if error;
% param_override.cluster.type = 'torque';
% param_override.cluster.type = 'matlab';
param_override.cluster.type = 'debug';
% param_override.cluster.rerun_only = true;
% param_override.cluster.desired_time_per_job  = 240*60;
% param_override.cluster.cpu_time_mult  = 2;
% param_override.cluster.mem_mult  = 2;

%% Automated Section
% =====================================================================

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

% Process each of the segments
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  param.day_seg
  
  echogram_fn_dir = ct_filename_out(param,param.get_echogram_stats.data_type);
  
  % Load frames associated with this segment
  frames_fn = ct_filename_support(param,'','frames');
  load(frames_fn);
  
  % Load layer information
  layer_params = [];
  layer_params.name = 'surface';
  layer_params.source = 'layerdata';
  layer_params.existence_check = false;
  
  surf = opsLoadLayers(param,layer_params);

  noise_bins = param.get_echogram_stats.noise_bins;
  signal_bins = param.get_echogram_stats.signal_bins;
  
  % dt_frm
  % bins_frm: time in units of dt_frm
  gps_time = [];
  bins_frm = [];
  sums_frm = [];
  counts_frm = [];
  signal_max_vals = []
  signal_mean_vals = [];
  noise_max_vals = []
  noise_mean_vals = [];
  for frm = 1:10%length(frames.frame_idxs)
    if param.get_echogram_stats.echogram_img == 0
      echogram_fn = fullfile(echogram_fn_dir,sprintf('Data_%s_%03d.mat',param.day_seg,frm));
    else
      echogram_fn = fullfile(echogram_fn_dir, ...
        sprintf('Data_img_%02d_%s_%03d.mat',param.get_echogram_stats.echogram_img,param.day_seg,frm));
    end
    mdata = load_L1B(echogram_fn);
    if isempty(mdata.Time)
      continue;
    end
    dt_frm = mdata.Time(2)-mdata.Time(1);
    if frm == 1
      dt = dt_frm;
    elseif abs((dt_frm-dt)/dt_frm) > 1e-6
      error('dt has changed from %.14g to %14g.', dt, dt_frm);
    end
    Nt = size(mdata.Data,1);
    Nx = size(mdata.Data,2);
    
    % Surface
    surf_bins = interp1(surf.gps_time,surf.twtt,mdata.GPS_time);
    surf_bins = interp_finite(surf_bins,NaN);
    surf_bins = interp1(mdata.Time,1:Nt,surf_bins);
    
    % Create signal mask
    mask = false(Nt,Nx);
    for rline = 1:Nx
      if ~isnan(surf_bins(rline))
        mask(max(1,surf_bins(rline)+signal_bins(1)) : min(Nt,surf_bins(rline)+signal_bins(end)),rline) = true;
      end      
    end
    data = lp(mdata.Data);
    data(~isfinite(data)) = NaN;
    data(mask) = NaN;
    
    % Mean value results
    bins_frm = mdata.Time/dt_frm;
    sums_frm = nansum(data,2); % Take the mean of all valid samples
    counts_frm = sum(~isnan(data),2); % How many samples were used to calculate the mean
    
    if isempty(bins)
      bins = bins_frm;
      sums = zeros(size(sums_frm));
      counts = zeros(size(counts_frm));
    end
    if bins_frm(1) < bins(1)
      % Add new bins_frm to beginning of bins
      bins = bins_frm(1) : bins(end);
      sums = [zeros(bins(1)-bins_frm(1),1); sums];
      counts = [zeros(bins(1)-bins_frm(1),1); counts];
    end
    if bins_frm(end) > bins(end)
      % Add new bins_frm to end of bins
      bins = bins(1) : bins_frm(end);
      sums = [sums; zeros(bins_frm(end)-bins(end),1)];
      counts = [counts; zeros(bins_frm(end)-bins(end),1)];
    end
    bin_start = bins_frm(1)-bins(1);
    sums(bin_start + (1:Nt)) = sums(bin_start + (1:Nt)) + sums_frm;
    counts(bin_start + (1:Nt)) = counts(bin_start + (1:Nt)) + counts_frm;

    % Concatenated results
    % 1. gps_time
    gps_time(end+(1:Nx)) = mdata.GPS_time;
    
    % 2. signal
    data = lp(mdata.Data);
    data(~isfinite(data)) = NaN;
    data(~mask) = NaN;
    cur_rline = length(signal_max_vals);
    signal_max_vals(end+(1:Nx)) = zeros(1,Nx);
    signal_mean_vals(end+(1:Nx)) = zeros(1,Nx);
    for rline = 1:Nx
      signal_max_vals(cur_rline+rline) = nanmax(data(:,rline),[],1);
      signal_mean_vals(cur_rline+rline) = nanmean(data(:,rline),1);
    end
    
    % 3. signal
    % Create noise mask
    mask = false(Nt,Nx);
    for rline = 1:Nx
      if ~isnan(surf_bins(rline))
        mask(max(1,surf_bins(rline)+noise_bins(1)) : min(Nt,surf_bins(rline)+noise_bins(end)),rline) = true;
      end      
    end
    data = lp(mdata.Data);
    data(~isfinite(data)) = NaN;
    data(~mask) = NaN;
    cur_rline = length(noise_max_vals);
    noise_max_vals(end+(1:Nx)) = zeros(1,Nx);
    noise_mean_vals(end+(1:Nx)) = zeros(1,Nx);
    for rline = 1:Nx
      noise_max_vals(cur_rline+rline) = nanmax(data(:,rline),[],1);
      noise_mean_vals(cur_rline+rline) = nanmean(data(:,rline),1);
    end
  end
  
end

