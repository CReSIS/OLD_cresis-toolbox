function data_pulse_compress

%% Pre-pulse compression filter
% ===================================================================

%% Burst RFI removal
% ===================================================================

%% Load coherent noise
% ===================================================================
if strcmpi(param.proc.coh_noise_method,'analysis')
  cdf_fn_dir = fileparts(ct_filename_out(param,'analysis', ''));
  cdf_fn = fullfile(cdf_fn_dir,sprintf('coh_noise_simp_%s.nc', param.day_seg));
  
  finfo = ncinfo(cdf_fn);
  % Determine number of records and set recs(1) to this
  Nt = finfo.Variables(find(strcmp('coh_aveI',{finfo.Variables.Name}))).Size(2);
  
  noise = [];
  noise.gps_time = ncread(cdf_fn,'gps_time');
  recs = find(noise.gps_time > records.gps_time(1) - 100 & noise.gps_time < records.gps_time(end) + 100);
  noise.gps_time = noise.gps_time(recs);
  
  noise.coh_ave = ncread(cdf_fn,'coh_aveI',[recs(1) 1],[recs(end)-recs(1)+1 Nt]) ...
    + j*ncread(cdf_fn,'coh_aveQ',[recs(1) 1],[recs(end)-recs(1)+1 Nt]);
  
  noise.coh_ave = - interp1(reshape(noise.gps_time,[numel(noise.gps_time) 1]),noise.coh_ave,records.gps_time,'linear','extrap').';
elseif strcmpi(param.get_heights.coh_noise_method,'estimate')
  % Apply coherent noise methods that require estimates derived now
  
  if param.get_heights.coh_noise_arg.DC_remove_en
    for wf_adc_idx = 1:size(g_data,3)
      g_data(:,:,wf_adc_idx) = bsxfun(@minus, g_data(:,:,wf_adc_idx), ...
        mean(g_data(:,:,wf_adc_idx),2));
    end
  end
  
  if length(param.get_heights.coh_noise_arg.B_coh_noise) > 1
    if length(param.get_heights.coh_noise_arg.A_coh_noise) > 1
      % Use filtfilt
      for wf_adc_idx = 1:size(g_data,3)
        g_data(:,:,wf_adc_idx) = single(filtfilt(param.get_heights.coh_noise_arg.B_coh_noise, ...
          param.get_heights.coh_noise_arg.A_coh_noise, double(g_data(:,:,wf_adc_idx).'))).';
      end
    else
      % Use fir_dec (no feedback)
      for wf_adc_idx = 1:size(g_data,3)
        g_data(:,:,wf_adc_idx) = fir_dec(g_data(:,:,wf_adc_idx),param.get_heights.coh_noise_arg.B_coh_noise,1);
      end
    end
  end
  
end




%% DDC, filter, fast-time-decimation
% ===================================================================

%% Pulse compression
% ===================================================================
if strcmpi(param.proc.pulse_compress,'matched')
elseif strcmpi(param.proc.pulse_compress,'deramp')
elseif strcmpi(param.proc.pulse_compress,'stepped')
end

%% Oversampling
% ===================================================================
if ~isequal(param.proc.ft_oversample,[1 1])
  for img = 1:length(data)
    for wf_adc_idx = 1:size(data,3)
      data = resample(data, param.proc.ft_oversample(1), param.proc.ft_oversample(2));
    end
  end
end

end


