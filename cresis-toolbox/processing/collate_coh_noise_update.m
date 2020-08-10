function collate_coh_noise_update(param,param_override)
% collate_coh_noise_update(param,param_override)
%
% Update collate_coh_noise files
%
% Example:
%  See run_collate_coh_noise_update for how to run.
%
% Authors: John Paden

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Update files if needed
% =====================================================================
if cluster_job_check()
  error('collate_coh_noise_update may not be called from cluster_job (gRadar.cluster.is_cluster_job is currently set to true). To remove this error, run collate_coh_noise_update on:\n  %s\n  %s\n  %s', param.radar_name, param.season_name, param.day_seg);
end

% Loop through each waveform for this segment
for wf = 1:length(param.radar.wfs)
  % Make sure param.radar.wfs.rx_paths exists
  if ~isfield(param.radar.wfs,'rx_paths')
    param.radar.wfs.rx_paths = 1;
  end
  % Loop through each adc for this segment
  for adc = 1:length(param.radar.wfs(wf).rx_paths)
    if ~isfinite(param.radar.wfs(wf).rx_paths(adc))
      % NaN in rx_paths means the adc does not exist
      continue;
    end
    
    % Create filename
    noise_fn_dir = fileparts(ct_filename_out(param,param.radar.wfs(wf).coh_noise_arg.fn, ''));
    noise_fn = fullfile(noise_fn_dir,sprintf('coh_noise_simp_%s_wf_%d_adc_%d.mat', param.day_seg, wf, adc));

    % Check for existence
    if ~exist(noise_fn,'file')
      warning('File does not exist so skipping: %s', noise_fn);
      continue;
    end
    
    % Load file
    fprintf('  Load coh_noise: %s (%s)\n', noise_fn, datestr(now));
    noise = load(noise_fn);
    
    % Check to see if missing any version 1 fields
    update_collate_coh_noise_file = false;
    if ~isfield(noise,'param_collate_coh_noise') || isempty(noise.param_collate_coh_noise)
      update_collate_coh_noise_file = true;
      warning('Missing field. Updated file will be saved.');
      noise.param_collate_coh_noise = noise.param_collate;
      noise = rmfield(noise,'param_collate');
    end
    if ~isfield(noise,'param_records') || isempty(noise.param_records)
      update_collate_coh_noise_file = true;
      warning('Missing field. Updated file will be saved.');
      noise.param_records = param;
    end
    if ~isfield(noise.param_collate_coh_noise.collate_coh_noise,'method') || isempty(noise.param_collate_coh_noise.collate_coh_noise.method)
      update_collate_coh_noise_file = true;
      warning('Missing field. Updated file will be saved.');
      noise.param_collate_coh_noise.collate_coh_noise.method = 'dft';
    end
    if isfield(noise,'coh_noise') && isfield(noise,'coh_noise_gps_time')
      update_collate_coh_noise_file = true;
      warning('Missing field. Updated file will be saved.');
      noise.firdec_gps_time = noise.coh_noise_gps_time;
      noise.firdec_noise    = noise.coh_noise;
      noise = rmfield(noise,{'coh_noise','coh_noise_gps_time'});
    end
    if ~isfield(noise.param_analysis.radar.wfs(wf),'system_dB')
      update_collate_coh_noise_file = true;
      warning('Missing field. Updated file will be saved.');
      % Update each entry with system_dB default value of 0
      for tmp_wf = 1:length(param.radar.wfs)
        noise.param_analysis.radar.wfs(tmp_wf).system_dB = 0;
      end
    end
      
    % Save result if missing any version 1 fields
    if update_collate_coh_noise_file
      ct_save(noise_fn,'-struct','noise')
    end
  end
end
