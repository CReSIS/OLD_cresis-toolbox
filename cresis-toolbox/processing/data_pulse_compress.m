function [data] = data_pulse_compress(param,hdr,wfs,data)
% [data] = data_pulse_compress(param,hdr,wfs,data)

[output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);

for img = 1:length(param.load.imgs)
  for wf_adc_idx = 1:size(data{img},3)
    wf = param.load.imgs{img}(wf_adc_idx,1);
    adc = param.load.imgs{img}(wf_adc_idx,2);
    
    %% Pre-pulse compression filter
    % ===================================================================
    
    %% Burst RFI removal
    % ===================================================================
    
    %% Load coherent noise
    % ===================================================================
    if strcmpi(wfs(wf).coh_noise_method,'analysis')
      cdf_fn_dir = fileparts(ct_filename_out(param,wfs(wf).coh_noise_arg.fn, ''));
      cdf_fn = fullfile(cdf_fn_dir,sprintf('coh_noise_simp_%s_wf_%d_adc_%d.nc', param.day_seg, wf, adc));
      
      finfo = ncinfo(cdf_fn);
      % Determine number of records and set recs(1) to this
      Nt = finfo.Variables(find(strcmp('coh_aveI',{finfo.Variables.Name}))).Size(2);
      
      noise = [];
      noise.gps_time = ncread(cdf_fn,'gps_time');
      recs = find(noise.gps_time > records.gps_time(1) - 100 & noise.gps_time < records.gps_time(end) + 100);
      noise.gps_time = noise.gps_time(recs);
      
      noise.coh_ave = ncread(cdf_fn,'coh_aveI',[recs(1) 1],[recs(end)-recs(1)+1 Nt]) ...
        + 1i*ncread(cdf_fn,'coh_aveQ',[recs(1) 1],[recs(end)-recs(1)+1 Nt]);
      
      data{img}(1:wfs(wf).Nt_raw,:,wf_adc_idx) = data{img}(1:wfs(wf).Nt_raw,:,wf_adc_idx) ...
        -interp1(reshape(noise.gps_time,[numel(noise.gps_time) 1]),noise.coh_ave,records.gps_time,'linear','extrap').';
      
    elseif strcmpi(wfs(wf).coh_noise_method,'estimated')
      % Apply coherent noise methods that require estimates derived now
      
      if wfs(wf).coh_noise_arg.DC_remove_en
        data{img}(1:wfs(wf).Nt_raw,:,wf_adc_idx) = bsxfun(@minus, data{img}(1:wfs(wf).Nt_raw,:,wf_adc_idx), ...
          mean(data{img}(1:wfs(wf).Nt_raw,:,wf_adc_idx),2));
      end
      
      if length(wfs(wf).coh_noise_arg.B_coh_noise) > 1
        if length(wfs(wf).coh_noise_arg.A_coh_noise) > 1
          % Use filtfilt (feedback)
          data{img}(1:wfs(wf).Nt_raw,:,wf_adc_idx) = single(filtfilt(wfs(wf).coh_noise_arg.B_coh_noise, ...
            wfs(wf).coh_noise_arg.A_coh_noise, double(data{img}(1:wfs(wf).Nt_raw,:,wf_adc_idx).'))).';
        else
          % Use fir_dec (no feedback)
          data{img}(1:wfs(wf).Nt_raw,:,wf_adc_idx) = fir_dec(data{img}(1:wfs(wf).Nt_raw,:,wf_adc_idx),wfs(wf).coh_noise_arg.B_coh_noise,1);
        end
      end
      
    end
    
    %% Pulse compression
    % ===================================================================
    if param.load.pulse_comp == 1
      
      if strcmpi(radar_type,'pulsed')
        if ~param.load.raw_data
          if param.load.pulse_comp
            % Digital down conversion
            data{img}(1:wfs(wf).Nt_raw,:,wf_adc_idx) = bsxfun(@times,data{img}(1:wfs(wf).Nt_raw,:,wf_adc_idx), ...
              exp(-1i*2*pi*(wfs(wf).fc-wfs(wf).DDC_freq)*wfs(wf).time_raw));
            
            % Pulse compression
            %   Apply matched filter and transform back to time domain
            tmp_data = circshift(ifft(bsxfun(@times,fft(data{img}(1:wfs(wf).Nt_raw,:,wf_adc_idx), wfs(wf).Nt_pc),wfs(wf).ref{adc})),wfs(wf).pad_length,1);
            
            % Decimation
            data{img}(1:wfs(wf).Nt,:,wf_adc_idx) = single(resample(double(tmp_data), wfs(wf).ft_dec(1), wfs(wf).ft_dec(2)));
          end
        end
        
      elseif strcmpi(radar_type,'deramp')
        
        wfs(wf).BW_window = [2.7e9 8e9];
        
        % Minimum step size in BW using only FFTs to constrain BW
        wfs(wf).df = wfs(wf).chirp_rate * 2^(1+wfs(wf).DDC_mode_max) / wfs(wf).fs_raw;
        
        % Check that number of samples is an integer
        wfs(wf).Nt = abs(diff(wfs(wf).BW_window)) / wfs(wf).df;
        if mod(wfs(wf).Nt,1)/wfs(wf).Nt > 1e3*eps
          error('BW_window bandwidth %g should be chosen so that number of samples (%g) for largest DDC_mod is an integer. It must be a multiple of %.16g Hz.', diff(wfs(wf).BW_window), wfs(wf).Nt, wfs(wf).df);
        end
        wfs(wf).Nt = round(wfs(wf).Nt) + 1;
        % Recalculate the end bandwidth point to reduce rounding errors
        wfs(wf).BW_window(2) = wfs(wf).BW_window(1) + (wfs(wf).Nt-1)*wfs(wf).df;
        
        for rec = 1:size(data{img},2)
          % Check to see if axes has changed since last record
          if rec == 1 ...
              || hdr.DDC_mode{img}(rec) ~= hdr.DDC_mode{img}(rec-1) ...
              || hdr.DDC_freq{img}(rec) ~= hdr.DDC_freq{img}(rec-1) ...
              || hdr.nyquist_zone{img}(rec) ~= hdr.nyquist_zone{img}(rec-1)
            
            if hdr.DDC_mode{img}(rec)
              fs_raw = wfs(wf).fs_raw ./ 2^(1+hdr.DDC_mode{img}(rec));
              Nt = wfs(wf).Nt * 2.^(-hdr.DDC_mode{img}(rec)+wfs(wf).DDC_mode_max);
            else
              fs_raw = wfs(wf).fs_raw;
              Nt = wfs(wf).Nt * 2.^(1+wfs(wf).DDC_mode_max);
            end
            
            t0 = hdr.t0{img}(rec) + wfs(wf).Tadc_adjust;
            Nt_raw = hdr.Nt{img}(rec);
            dt_raw = 1/fs_raw;
            time_raw = (t0:dt_raw:t0+dt_raw*(Nt_raw-1)).';
            
            df_raw = 1/(dt_raw*Nt_raw);
            
            DDC_freq = hdr.DDC_freq{img}(rec);
            
            nz = double(hdr.nyquist_zone{img}(rec));
            
            f_nz0 = wfs(wf).fs_raw * floor(nz/2);
            
            freq_raw =  f_nz0 + mod(DDC_freq + df_raw*ifftshift(-floor(Nt_raw/2):floor((Nt_raw-1)/2)).', wfs(wf).fs_raw);
            freq_raw_valid = freq_raw;
            
            conjugate_bins = ~(freq_raw_valid >= nz*wfs(wf).fs_raw/2 ...
              & freq_raw_valid <= (1+nz)*wfs(wf).fs_raw/2);
            
            if mod(nz,2)
              freq_raw_valid(conjugate_bins) = nz*wfs(wf).fs_raw-freq_raw_valid(conjugate_bins);
            else
              freq_raw_valid(conjugate_bins) = (nz+1)*wfs(wf).fs_raw - freq_raw_valid(conjugate_bins);
            end
            
            freq_raw_valid = df_raw*round(freq_raw_valid/df_raw);
            
            [~,unique_idxs] = unique(freq_raw_valid);
            
            freq_raw_unique = freq_raw_valid(unique_idxs);
            conjugate_unique = conjugate_bins(unique_idxs);
            % unique_idxs(~valid_bins);
            
            if 0
              figure(1);
              clf;
              plot(freq_raw_valid,'.')
              hold on
              freq_raw_valid(conjugate_bins) = NaN;
              plot(freq_raw_valid,'r.')
              grid on
              
              figure(2);
              clf;
              plot(freq_raw_unique,'.')
              hold on
              tmp = freq_raw_unique;
              tmp(~conjugate_unique) = NaN;
              plot(tmp,'r.')
              grid on
            end
            
            % Create window
            H_Nt = wfs(wf).ft_wind(Nt);
            
            
            %             wfs(wf).fs = 1/dt;
            %             wfs(wf).Nt = length(time);
            %             wfs(wf).freq = wfs(wf).fc + ifftshift( -floor(wfs(wf).Nt/2)*df : df : floor((wfs(wf).Nt-1)/2)*df ).';
            
            if 0
              figure(1); clf;
              plot(time);
            end
            
          end
          
          % Check to see if axes has changed since last record
          if rec == 1 ...
              || hdr.t_ref{img}(rec) ~= hdr.t_ref{img}(rec-1) ...
              
            % Convert IF frequency to time delay
            time = freq_raw_unique/wfs(wf).chirp_rate - hdr.t_ref{img}(rec);
            
            % Ensure that start time is a multiple of dt
            dt = time(2)-time(1);
            time_correction = dt - mod(time(1),dt);
            time = time + time_correction;
            time_correction_IF = -time_correction*wfs(wf).chirp_rate;
          end
          
          % Convert raw time into instantaneous frequency for the surface bin
          f_rf = wfs(wf).f0 + wfs(wf).chirp_rate*(time_raw - hdr.surface(rec));
          if wfs(wf).BW_window(2) > max(f_rf)
            error('BW_window (%g) is more than maximum measured RF frequency (%g) with surface twtt %g.', wfs(wf).BW_window(2), max(f_rf), hdr.surface(rec))
          end
          if wfs(wf).BW_window(1) < min(f_rf)
            error('BW_window (%g) is less than minimum measured RF frequency (%g) with surface twtt %g.', wfs(wf).BW_window(1), min(f_rf), hdr.surface(rec))
          end
          
          % Create the window for the particular range line
          window_start_idx = find(f_rf >= wfs(wf).BW_window(1),1);
          H = zeros(size(f_rf));
          H(window_start_idx : window_start_idx+Nt-1) = H_Nt;
          
          % Modulate the raw data to adjust the start time to always be a
          % multiple of wfs(wf).dt
          hdr.t0{img}(rec) = time(1);
          data{img}(1:hdr.Nt{img}(rec),rec,wf_adc_idx) = data{img}(1:hdr.Nt{img}(rec),rec,wf_adc_idx) .* exp(1i*2*pi*time_correction_IF*time_raw);
          
          % Window and Pulse compress
          if 0
            clf;
            plot(lp(data{img}(1:hdr.Nt{img}(rec),rec,wf_adc_idx)));
            hold on;
            plot(lp(data{img}(1:hdr.Nt{img}(rec),rec,wf_adc_idx) .* H));
          end
          tmp = fft(data{img}(1:hdr.Nt{img}(rec),rec,wf_adc_idx) .* H);
          % Reorder result in case it is wrapped
          tmp = tmp(unique_idxs);
          % Some of the frequency bins are conjugated versions of the
          % signal
          tmp(conjugate_unique) = conj(tmp(conjugate_unique));
          
          % Complex baseband data
          tmp = ifft(ifftshift(fft(conj(tmp))));
          
          % Deskew of the residual video phase
          deskew = exp(-1i*(pi*wfs(wf).chirp_rate*(time.^2 - hdr.t_ref{img}(rec).^2)));
          tmp = tmp.*deskew;
          
          % Update the data matrix with the pulse compressed waveform
          hdr.Nt{img}(rec) = length(tmp);
          data{img}(1:hdr.Nt{img}(rec),rec,wf_adc_idx) = tmp;

          if 0
            fs = 125e9;
            Tpd = 240e-6;
            f0 = 2e9;
            f1 = 18e9;
            fs_IF = 250e6;
            
            dt = 1/fs;
            Nt = fs*Tpd;
            alpha = (f1-f0)/Tpd;
            fc = (f0+f1)/2;
            time = dt*(0:Nt-1).';
            s = cos(2*pi*f0*time + pi*alpha*time.^2);
            BW = f1-f0;
            sigma_t = 1/BW;
            delta_ts = linspace(0,sigma_t,11);
            t0 = 1e-6;
            for delta_t = delta_ts
              td = t0+delta_t;
              sd = fft(s.*cos(2*pi*f0*(time-td) + pi*alpha*(time-td).^2));
              df = 1/(dt*Nt);
              freq_raw = df * (0:Nt-1).';
              Mt = fs/fs_IF;
              sd = sd(freq_raw < fs_IF);
%               sd = ifft(conj(ifftshift(fft(sd))));
              time_IF = freq_raw(freq_raw < fs_IF)/alpha;
              %             sd = resample(sd,Mt,1);
              clf;
              plot(time_IF, lp(sd));
              angle(max(sd))
              angle(exp(1i*(2*pi*(f0+alpha*time)*(-td) - pi*alpha*td.^2 - pi*alpha*t0.^2)))
              %2*pi*f0*delta_t
              
              xlim([t0 + [-2 2]*sigma_t])
              td
              pause
            end
            
          end
          
          if 0
            clf;
            plot(time, lp(data{img}(1:hdr.Nt{img}(rec),rec,wf_adc_idx)));
          end
        end
        
        % Create a matrix of data with constant time rows, fill invalid samples with NaN
        tmp = data{img}(1:hdr.Nt{img}(rec),:,wf_adc_idx);
        tmp = interpft(tmp,10*size(tmp,1));
        imagesc(lp(tmp))
        
        [dd,dd2]=max(tmp);
        physical_constants;
        dd = dd.*exp(1i*2*pi*hdr.records{1}.elev / (c/2) * mean(wfs(wf).BW_window));
        plot(dd2 - mean(dd2))
        plot(lp(dd))
        plot(angle(dd))
        
        keyboard
        
      elseif strcmpi(radar_type,'stepped')
        
      end
    end
    
    % %% Oversampling <-- Should be done after all processing is completed
    % % ===================================================================
    % if ~isequal(param.load.ft_oversample,[1 1])
    %   for img = 1:length(data)
    %     for wf_adc_idx = 1:size(data,3)
    %       data = resample(data, param.load.ft_oversample(1), param.load.ft_oversample(2));
    %     end
    %   end
    % end
    
  end
  data{img} = data{img}(1:wfs(wf).Nt,:,:);
end
