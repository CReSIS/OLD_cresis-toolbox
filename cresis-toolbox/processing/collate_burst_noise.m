% function collate_burst_noise(param,param_override)
% collate_burst_noise(param,param_override)
%
% Collects analysis.m results from burst noise tracking (burst_noise
% command) and creates files for removing the burst noise during data
% loading.
% Loads all the burst_noise_* files and creates burst_noise_simp_* files that
% are pre-filtered for speed and saves as netcdf so that subsets of the
% files can be loaded efficiently.
%
% Example:
%  See run_collate_burst_noise for how to run.
%
% Authors: John Paden

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input checks
% =====================================================================

if ~isfield(param,'collate_burst_noise') || isempty(param.collate_burst_noise)
  param.collate_burst_noise = [];
end

% param.collate_burst_noise.imgs: Check this first since other input checks
%   are dependent on its value.
if ~isfield(param.analysis,'imgs') || isempty(param.analysis.imgs)
  param.analysis.imgs = {[1 1]};
end
if ~isfield(param.collate_burst_noise,'imgs') || isempty(param.collate_burst_noise.imgs)
  param.collate_burst_noise.imgs = 1:length(param.analysis.imgs);
end

if ~isfield(param.collate_burst_noise,'cmd_idx') || isempty(param.collate_burst_noise.cmd_idx)
  param.collate_burst_noise.cmd_idx = 1;
end

if ~isfield(param.collate_burst_noise,'debug_plots') || isempty(param.collate_burst_noise.debug_plots)
  param.collate_burst_noise.debug_plots = {'bn_plot'};
end
enable_visible_plot = any(strcmp('visible',param.collate_burst_noise.debug_plots));
enable_bn_plot = any(strcmp('bn_plot',param.collate_burst_noise.debug_plots));
if ~isempty(param.collate_burst_noise.debug_plots)
  h_fig = get_figures(2,enable_visible_plot);
end

if ~isfield(param.(mfilename),'debug_out_dir') || isempty(param.(mfilename).debug_out_dir)
  param.(mfilename).debug_out_dir = mfilename;
end
debug_out_dir = param.(mfilename).debug_out_dir;

if ~isfield(param.collate_burst_noise,'in_path') || isempty(param.collate_burst_noise.in_path)
  param.collate_burst_noise.in_path = 'analysis';
end

if ~isfield(param.collate_burst_noise,'max_bad_waveforms') || isempty(param.collate_burst_noise.max_bad_waveforms)
  param.collate_burst_noise.max_bad_waveforms = 1;
end

if ~isfield(param.collate_burst_noise,'out_path') || isempty(param.collate_burst_noise.out_path)
  param.collate_burst_noise.out_path = param.collate_burst_noise.in_path;
end

if ~isfield(param.collate_burst_noise,'wf_adcs') || isempty(param.collate_burst_noise.wf_adcs)
  param.collate_burst_noise.wf_adcs = [];
end
if ~isempty(param.collate_burst_noise.wf_adcs) && ~iscell(param.collate_burst_noise.wf_adcs)
  wf_adcs = param.collate_burst_noise.wf_adcs;
  param.collate_burst_noise.wf_adcs = {};
  for img = 1:max(param.collate_burst_noise.imgs)
    param.collate_burst_noise.wf_adcs{img} = wf_adcs;
  end
end

records = records_load(param);
new_bit_mask = zeros(size(records.bit_mask));

for img = param.collate_burst_noise.imgs
  
  if isempty(param.collate_burst_noise.wf_adcs)
    wf_adcs = 1:size(param.analysis.imgs{img},1);
  else
    wf_adcs = param.collate_burst_noise.wf_adcs{img};
  end
  for wf_adc = wf_adcs
    wf = param.analysis.imgs{img}(wf_adc,1);
    adc = param.analysis.imgs{img}(wf_adc,2);
    
    %% Load the burst noise file
    % ===================================================================
    fn_dir = fileparts(ct_filename_out(param,param.collate_burst_noise.in_path));
    fn = fullfile(fn_dir,sprintf('burst_noise_%s_wf_%d_adc_%d.mat', param.day_seg, wf, adc));
    fprintf('Loading %s (%s)\n', fn, datestr(now));
    noise = load(fn);
    
    cmd = noise.param_analysis.analysis.cmd{param.collate_burst_noise.cmd_idx};
    
    bad_recs_unique = unique(noise.bad_recs);
    
    % Determine the boards for this wf-adc pair
    param.load.imgs = {[wf adc]};
    [wfs,states] = data_load_wfs(param, records);
    % Combine all boards: Some wf-adc pairs result in multiple waveforms
    % being loaded (e.g. for IQ on transmit or separate IQ channels or
    % zero-pi mod sequences that are stored separately)
    boards = cell2mat({states.board_idx});
    
    % Set bit 2 to true for the bad records
    new_bit_mask(boards,bad_recs_unique) = bitor(4,new_bit_mask(boards,bad_recs_unique));
    
    %% Plot
    % =====================================================================
    if enable_bn_plot
      clf(h_fig(1));
      set(h_fig(1), 'name', 'burst_noise rec-bin');
      h_axes(1) = axes('parent',h_fig(1));
      plot(noise.bad_recs, noise.bad_bins, 'x', 'parent', h_axes(1));
      title(h_axes(1), sprintf('%s wf %d adc %d',regexprep(param.day_seg,'_','\\_'), wf, adc));
      xlabel(h_axes(1), 'Record');
      ylabel(h_axes(1), 'Range bin');
      
      fig_fn = [ct_filename_ct_tmp(param,'',debug_out_dir,sprintf('burst_rec_bin_wf_%02d_adc_%02d',wf,adc)) '.jpg'];
      fprintf('Saving %s\n', fig_fn);
      fig_fn_dir = fileparts(fig_fn);
      if ~exist(fig_fn_dir,'dir')
        mkdir(fig_fn_dir);
      end
      ct_saveas(h_fig(1),fig_fn);
      fig_fn(end-2:end) = 'fig';
      fprintf('Saving %s\n', fig_fn);
      ct_saveas(h_fig(1),fig_fn);
      
      clf(h_fig(2));
      set(h_fig(2), 'name', 'burst_noise waveforms');
      h_axes(2) = axes('parent',h_fig(2));
      for block_idx = 1:length(noise.bad_waveforms)
        plot(lp(noise.bad_waveforms{block_idx}(:,1:min(param.collate_burst_noise.max_bad_waveforms,end))), 'parent', h_axes(2));
        hold(h_axes(2),'on');
      end
      title(h_axes(2), sprintf('%s wf %d adc %d',regexprep(param.day_seg,'_','\\_'), wf, adc));
      xlabel(h_axes(2), 'Range bin');
      ylabel(h_axes(2), 'Relative power (dB)');
      
      fig_fn = [ct_filename_ct_tmp(param,'',debug_out_dir,sprintf('burst_waveform_wf_%02d_adc_%02d',wf,adc)) '.jpg'];
      fprintf('Saving %s\n', fig_fn);
      fig_fn_dir = fileparts(fig_fn);
      if ~exist(fig_fn_dir,'dir')
        mkdir(fig_fn_dir);
      end
      ct_saveas(h_fig(2),fig_fn);
      fig_fn(end-2:end) = 'fig';
      fprintf('Saving %s\n', fig_fn);
      ct_saveas(h_fig(2),fig_fn);
    end
 
    if enable_visible_plot && ~isempty(bad_recs_unique)
      % Bring plots to front
      for h_fig_idx = 1:length(h_fig)
        figure(h_fig(h_fig_idx));
      end
      % Enter debug mode
      keyboard
    end
    
  end
end

%% Update records file
records.bit_mask = records.bit_mask - bitand(4,records.bit_mask) + uint8(new_bit_mask);
records_fn = ct_filename_support(param,'','records');
ct_save(records_fn,'-append','-struct','records','bit_mask');

if ~enable_visible_plot
  try
    delete(h_fig);
  end
end
