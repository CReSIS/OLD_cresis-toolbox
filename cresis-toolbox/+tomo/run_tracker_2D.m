%% Script run_tracker_2D
%
% Runs the selected 2D tracking algorithms
%  on the desired dataset
%
% Authors: Victor Berger
%
% See also: viterbi_tracker_2D.m, mcmc_tracker_2D.m

dbstack_info = dbstack;
fprintf('=====================================================================\n');
fprintf('%s (%s)\n', dbstack_info(1).name, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

%% General User Settings
% Tracking algorithms: 'viterbi', 'mcmc', 'lsm'
% algorithms = {'viterbi', 'mcmc', 'lsm'};
algorithms = {'viterbi'};

% params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'','post');
% params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20140325_05');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20140325_06');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20140325_07');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20140401_03');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20140506_01');


params = read_param_xls(ct_filename_param('rds_param_2018_Greenland_P3.xls'),'','post');
params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20180404_02');
params = ct_set_params(params,'cmd.generic',1,'day_seg','20180405_01');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20180406_01');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20180418_04');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20180418_05');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20180418_06');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20180419_01');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20180419_02');

params = ct_set_params(params,'cmd.frms',[]);

options.name       = 'mvdr';
options.debug      = false;
options.ops_write  = true;

options.surf_layer.name = 'surface';
options.surf_layer.source = 'layerData';

%% Ice mask options
if 1 % If using GeoTIFF file for ice mask
  options.binary_icemask = false;
  options.icemask_fn = 'greenland/IceMask/GimpIceMask_90m_v1.1.tif';
  
  % Useful for Antarctica seasons:
  % options.icemask_fn  = 'antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_icemask_grounded_and_shelves.tif';
  % options.icemask2_fn = 'antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_rockmask.tif';
  
  if ~isempty(options.icemask_fn)
    options.icemask_fn  = ct_filename_gis([],options.icemask_fn);
  end
  if isfield(options, 'icemask2_fn') && ~isempty(options.icemask2_fn)
    options.icemask2_fn = ct_filename_gis([],options.icemask2_fn);
  end
  
else % If using binary file for ice mask
  options.binary_icemask = true;
  param.radar_name       = 'rds';
  param.season_name      = '2014_Greenland_P3';
  options.icemask_fn     = ct_filename_gis(param,fullfile('canada','ice_mask','03_rgi50_ArcticCanadaNorth','03_rgi50_ArcticCanadaNorth.bin'));
end

%% Viterbi User Settings
options.viterbi.crossoverload  = true;
options.viterbi.layername      = 'viterbi_bottom';
options.viterbi.framecat       = false;

options.viterbi.layer_dest.name = 'bottom';
options.viterbi.layer_dest.source = 'layerdata';
% options.viterbi.layer_dest.layerdata_source = 'layerdata';

% options.viterbi.layer_dest.source = 'ops';
% options.viterbi.layer_dest.group = 'standard';
% options.viterbi.layer_dest.description = '';

options.viterbi.bottom_bin     = -1;
options.viterbi.egt_weight     = -1;
options.viterbi.mu_size        = 31;
options.viterbi.mu             = log10(exp(-(-(options.viterbi.mu_size-1)/2 : (options.viterbi.mu_size-1)/2).^4/1));
options.viterbi.mu_thr         = -30;
options.viterbi.mu(options.viterbi.mu < options.viterbi.mu_thr) = options.viterbi.mu_thr;
options.viterbi.mu             = options.viterbi.mu - mean(options.viterbi.mu);
options.viterbi.sigma          = sum(abs(options.viterbi.mu))/10*ones(1,options.viterbi.mu_size);
options.viterbi.smooth_var     = inf;
options.viterbi.repulsion      = 150000;
options.viterbi.smooth_weight  = 5;
options.viterbi.ice_bin_thr    = 10;
options.viterbi.CF.sensorydist = 200;
options.viterbi.CF.max_cost    = 50;
options.viterbi.CF.lambda      = 0.075;

%% MCMC User Settings
options.mcmc.alg           = 'HMM'; % 'MCMC' or 'HMM'
options.mcmc.lyrtop        = 'mcmc_top';
options.mcmc.lyrbot        = 'mcmc_bot';

options.mcmc.top_smooth    = 1000;
options.mcmc.bottom_smooth = 1000;
options.mcmc.top_peak      = 0.5;
options.mcmc.bottom_peak   = 0.5;
options.mcmc.repulsion     = 10;

%% LSM User Settings
options.lsm.lyrtop        = 'lsm_top';
options.lsm.lyrbot        = 'lsm_bot';

options.lsm.numOuterIter  = 175;
options.lsm.maxOuterIter  = 250;

%% Automated section
global gRadar;

for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  % Load frames file
  load(ct_filename_support(param,'','frames'));
  
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
  
  data_fn_dir = ct_filename_out(param, options.name, '');
  for frm = param.cmd.frms
    data_fn_name = sprintf('Data_%s_%03d.mat',param.day_seg,frm);
    data_fn      = fullfile(data_fn_dir, data_fn_name);
    try
      data_struct.(sprintf('data_%s_%03d',param.day_seg,frm)) = load(data_fn);
    catch ME
      fprintf('\nFailed to load file %s, skipping.\n', data_fn);
      continue;
    end
  end
end

if any(strcmp(algorithms, 'viterbi'))
  [viterbi_layers, surf_bins] = tomo.viterbi_tracker_2D(params, options, data_struct);
end

if any(strcmp(algorithms, 'mcmc'))
  mcmc_layers = tomo.mcmc_tracker_2D(params, options, data_struct);
end

if any(strcmp(algorithms, 'lsm'))
  lsm_layers = tomo.lsm_tracker_2D(params, options, data_struct);
end

% Display results overlaid on echogram
if 0
  frame = '20170322_05_001';
  figure; imagesc(lp(data_struct.(sprintf('data_%s', frame)).Data));
  colormap(1-gray(256)); hold on;
  if any(strcmp(algorithms, 'viterbi'))
    plot(surf_bins);
    plot(viterbi_layers.(sprintf('layer_%s', frame)).bot);
  end
  if any(strcmp(algorithms, 'mcmc'))
    plot(mcmc_layers.(sprintf('layer_%s', frame)).bot);
    if 0
      plot(mcmc_layers.(sprintf('layer_%s', frame)).top);
    end
  end
  if any(strcmp(algorithms, 'lsm'))
    plot(lsm_layers.(sprintf('layer_%s', frame)).bot);
    if 0
      plot(lsm_layers.(sprintf('layer_%s', frame)).top);
    end
  end
  
end