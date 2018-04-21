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
  
        if ~param.proc.raw_data
        % Apply channel compensation
        chan_equal = 10.^(param.radar.wfs(wf).chan_equal_dB(param.radar.wfs(wf).rx_paths(adc))/20) ...
          .* exp(1i*param.radar.wfs(wf).chan_equal_deg(param.radar.wfs(wf).rx_paths(adc))/180*pi);
        state(board+1).data{accum_idx} = state(board+1).data{accum_idx}/chan_equal;
        state(board+1).data{accum_idx} = state(board+1).data{accum_idx}/wfs(wf).adc_gains(adc);
        
        if param.load.pulse_comp
          % Apply blank (only should enable if sidelobe problems present)
          state(board+1).data{accum_idx}(wfs(wf).time_raw-param.radar.wfs(wf).Tsys(adc) <= param.surface(rec) + wfs(wf).blank) = 0;
          
          % Digital down conversion and decimation
          state(board+1).data{accum_idx} = state(board+1).data{accum_idx}.*exp(-1i*2*pi*(wfs(wf).fc-wfs(wf).DDC_freq)*wfs(wf).time_raw);
          state(board+1).data{accum_idx} = resample(double(state(board+1).data{accum_idx}), param.wfs(wf).ft_dec(1), param.wfs(wf).ft_dec(2));
          
          % Zero pad front: (the standard)
          state(board+1).data{accum_idx} = fft([zeros(wfs(wf).pad_length,1); state(board+1).data{accum_idx}]);
          % Zero pad end: (debug only)
          %state(board+1).data{accum_idx} = fft(state(board+1).data{accum_idx}, wfs(wf).Nt_pc);
          
          % Pulse compression
          %   Apply matched filter and transform back to time domain
          state(board+1).data{accum_idx} = ifft(state(board+1).data{accum_idx} .* wfs(wf).ref{adc});
          
        end
        
      
      
      
      % Store wf-adc pair in output matrix
      data{img}(:,out_idx,wf_adc_idx) = state(board+1).data{accum_idx} / param.proc.presums;
    end
  
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


