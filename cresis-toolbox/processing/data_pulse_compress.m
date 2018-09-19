function [hdr,data] = data_pulse_compress(param,hdr,wfs,data)
% [hdr,data] = data_pulse_compress(param,hdr,wfs,data)

if param.load.raw_data && param.load.pulse_comp
  error('Pulse compression (param.load.pulse_comp) cannot be enabled with raw data loading (param.load.raw_data).');
end

hdr.custom = [];

physical_constants;

[output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);

for img = 1:length(param.load.imgs)
  for wf_adc = 1:size(param.load.imgs{img},1)
    wf = param.load.imgs{img}(wf_adc,1);
    adc = param.load.imgs{img}(wf_adc,2);
    BW_window_max_warning_printed = false;
    BW_window_min_warning_printed = false;
    
    %% Pre-pulse compression filter
    % ===================================================================
    if strcmpi(radar_type,'deramp')
      if strcmpi(wfs(wf).prepulse_H.type,'filtfilt')
        data{img}(:,:,wf_adc) = single(filtfilt(wfs(wf).prepulse_H.B,1,double(data{img}(:,:,wf_adc))));
      end
    end
    
    %% Burst RFI removal
    % ===================================================================
    
    %% Remove coherent noise for pulsed radar
    % ===================================================================
    if strcmpi(radar_type,'pulsed')
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
        
        data{img}(1:wfs(wf).Nt_raw,:,wf_adc) = data{img}(1:wfs(wf).Nt_raw,:,wf_adc) ...
          -interp1(reshape(noise.gps_time,[numel(noise.gps_time) 1]),noise.coh_ave,records.gps_time,'linear','extrap').';
        
      elseif strcmpi(wfs(wf).coh_noise_method,'estimated')
        % Apply coherent noise methods that require estimates derived now
        
        if wfs(wf).coh_noise_arg.DC_remove_en
          data{img}(1:wfs(wf).Nt_raw,:,wf_adc) = bsxfun(@minus, data{img}(1:wfs(wf).Nt_raw,:,wf_adc), ...
            mean(data{img}(1:wfs(wf).Nt_raw,:,wf_adc),2));
        end
        
        if length(wfs(wf).coh_noise_arg.B_coh_noise) > 1
          if length(wfs(wf).coh_noise_arg.A_coh_noise) > 1
            % Use filtfilt (feedback)
            data{img}(1:wfs(wf).Nt_raw,:,wf_adc) = single(filtfilt(wfs(wf).coh_noise_arg.B_coh_noise, ...
              wfs(wf).coh_noise_arg.A_coh_noise, double(data{img}(1:wfs(wf).Nt_raw,:,wf_adc).'))).';
          else
            % Use fir_dec (no feedback)
            data{img}(1:wfs(wf).Nt_raw,:,wf_adc) = fir_dec(data{img}(1:wfs(wf).Nt_raw,:,wf_adc),wfs(wf).coh_noise_arg.B_coh_noise,1);
          end
        end
        
      end
    end
    
    %% Pulse compression
    % ===================================================================
    if param.load.pulse_comp == 1
      
      if strcmpi(radar_type,'pulsed')
        % Digital down conversion
        data{img}(1:wfs(wf).Nt_raw,:,wf_adc) = bsxfun(@times,data{img}(1:wfs(wf).Nt_raw,:,wf_adc), ...
          exp(-1i*2*pi*(wfs(wf).fc-wfs(wf).DDC_freq)*wfs(wf).time_raw));
        
        % Pulse compression
        %   Apply matched filter and transform back to time domain
        tmp_data = circshift(ifft(bsxfun(@times,fft(data{img}(1:wfs(wf).Nt_raw,:,wf_adc), wfs(wf).Nt_pc),wfs(wf).ref{adc})),wfs(wf).pad_length,1);
        
        % Decimation
        data{img}(1:wfs(wf).Nt,:,wf_adc) = single(resample(double(tmp_data), wfs(wf).ft_dec(1), wfs(wf).ft_dec(2)));
        
        if wf_adc == 1
          hdr.time{img} = wfs(wf).time;
          hdr.freq{img} = wfs(wf).freq;
        end
        
        
      elseif strcmpi(radar_type,'deramp')
        %% Pulse compress: deramp
        if 0
          % ENABLE_FOR_DEBUG
          % Create simulated data
          clear;
          img = 1;
          rec = 1;
          wf = 1;
          wf_adc = 1;
          hdr.DDC_dec{img}(rec) = 1;
          hdr.DDC_freq{img}(rec) = 0e6;
          hdr.nyquist_zone_signal{img}(rec) = 1;
          wfs(wf).fs_raw = 250e6;
          wfs(wf).BW_window = [2.7e9 17.5e9];
          wfs(wf).f0 = 2e9;
          wfs(wf).f1 = 18e9;
          wfs(wf).Tpd = 240e-6;
          wfs(wf).td_mean = 3e-6;
          BW = wfs(wf).f1-wfs(wf).f0;
          wfs(wf).chirp_rate = BW/wfs(wf).Tpd;
          hdr.surface(rec) = 3e-6;
          hdr.t_ref{img}(rec) = 0e-6;
          wfs(wf).ft_wind = @hanning;
          tguard = 1e-6;
          td_max = 4e-6;
          td_min = 1e-6;
          tref = hdr.t_ref{img}(rec);
          hdr.t0_raw{img}(rec) = -tguard + min(td_max,tref);
          wfs(wf).Tadc_adjust = 0;
          fs_nyquist = max(wfs(wf).f0,wfs(wf).f1)*5;
          Mt_oversample = ceil(fs_nyquist/wfs(wf).fs_raw);
          fs_rf = Mt_oversample*wfs(wf).fs_raw;
          Nt_rf = round((wfs(wf).Tpd + max(td_max,tref) - min(td_min,tref) + 2*tguard)*fs_rf/Mt_oversample)*Mt_oversample;
          hdr.Nt{img}(rec) = Nt_rf / Mt_oversample;
          
          t0 = hdr.t0_raw{img}(1);
          fs_raw_dec = wfs(wf).fs_raw ./ hdr.DDC_dec{img}(1);
          dt_raw = 1/fs_raw_dec;
          time_raw_no_trim = (t0:dt_raw:t0+dt_raw*(hdr.Nt{img}(rec)-1)).';
          if 0
            tds = hdr.surface(rec) + 1/BW/5*(0:11);
            hdr.surface = tds;
            %hdr.surface(:) = hdr.surface(1);% Enable or disable this line to simulate errors in surface estimate
          else
            time_raw_no_trim_transition = hdr.surface(rec)+(wfs(wf).BW_window(1) - wfs(wf).f0)/wfs(wf).chirp_rate;
            idx = find(time_raw_no_trim>time_raw_no_trim_transition,1);
            
            % Update surface to lie on a IF sample boundary
            hdr.surface = time_raw_no_trim(idx)-(wfs(wf).BW_window(1)-wfs(wf).f0)/wfs(wf).chirp_rate;
            
            tds = hdr.surface(1) + 1/wfs(wf).fs_raw*(0 + [-100-1/3 -30 -10 -4/3 -1/3 0 1/3 4/3 10 30 100+1/3]);
            %tds = hdr.surface(1) + 1/BW/5*[-2 -1 0 1 2];
            hdr.surface = tds;
            hdr.surface(:) = hdr.surface(1);% Enable or disable this line to simulate errors in surface estimate
          end
          Tpd = wfs(wf).Tpd;
          f0 = wfs(wf).f0;
          f1 = wfs(wf).f1;
          alpha = wfs(wf).chirp_rate;
          fs_rf = Mt_oversample*wfs(wf).fs_raw;
          time = hdr.t0_raw{img}(rec) + 1/fs_rf * (0:Nt_rf-1).';
          for rec = 1:length(tds)
            fprintf('Simulating %d of %d\n', rec, length(tds));
            td = tds(rec);
            
            f_rf = wfs(wf).f0 + wfs(wf).chirp_rate*(time_raw_no_trim - hdr.surface(rec));
            window_start_idx = find(f_rf >= wfs(wf).BW_window(1),1);
            fprintf('  window_start_idx: %d\n', window_start_idx);
            
            hdr.DDC_dec{img}(rec) = hdr.DDC_dec{img}(1);
            hdr.DDC_freq{img}(rec) = hdr.DDC_freq{img}(1);
            hdr.nyquist_zone_signal{img}(rec) = hdr.nyquist_zone_signal{img}(1);
            hdr.t_ref{img}(rec) = hdr.t_ref{img}(1);
            hdr.t0_raw{img}(rec) = hdr.t0_raw{img}(1);
            hdr.Nt{img}(rec) = hdr.Nt{img}(1);
            
            s = tukeywin_cont((time-Tpd/2-td)/Tpd,0) .* cos(2*pi*f0*(time-td) + pi*alpha*(time-td).^2);
            r = tukeywin_cont((time-Tpd/2-tref)/Tpd,0) .* cos(2*pi*f0*(time-tref) + pi*alpha*(time-tref).^2);
            s_if = s.*r; clear r s;
            s_if_theory = 0.5*tukeywin_cont((time-Tpd/2-td/2-tref/2)/(Tpd-abs(tref-td)),0) ...
              .* (cos(2*pi*f0*(time-td) + pi*alpha*(time-td).^2 - (2*pi*f0*(time-tref) + pi*alpha*(time-tref).^2)) ...
              + cos(2*pi*f0*(time-td) + pi*alpha*(time-td).^2 + (2*pi*f0*(time-tref) + pi*alpha*(time-tref).^2)));
            
            [Bfilt,Afilt] = butter(6, wfs(wf).fs_raw*hdr.nyquist_zone_signal{img}(rec) / (fs_rf/2));
            s_if = filtfilt(Bfilt,Afilt,s_if);
            s_if_theory = 0.5*tukeywin_cont((time-Tpd/2-td/2-tref/2)/(Tpd-abs(tref-td)),0) ...
              .* cos(2*pi*f0*(time-td) + pi*alpha*(time-td).^2 - (2*pi*f0*(time-tref) + pi*alpha*(time-tref).^2));
            
            s_if = s_if(1:Mt_oversample:end);
            
            data{img}(:,rec,wf_adc) = s_if;
          end
          store_data = data;
          store_Nt = hdr.Nt{img};
          
        elseif 0
          % ENABLE_FOR_DEBUG
          % Use previously generated simulated data
          data = store_data;
          hdr.Nt{img} = store_Nt;
        end
        
        freq_axes_changed = false;
        for rec = 1:size(data{img},2)
          
          % Check to see if axes has changed since last record
          if rec == 1 ...
              || hdr.DDC_dec{img}(rec) ~= hdr.DDC_dec{img}(rec-1) ...
              || hdr.DDC_freq{img}(rec) ~= hdr.DDC_freq{img}(rec-1) ...
              || hdr.nyquist_zone_signal{img}(rec) ~= hdr.nyquist_zone_signal{img}(rec-1)
            
            freq_axes_changed = true;
            
            %% Set output time sampling
            % The output time axes for every choice of DDC_dec must have
            % the same sample spacing. We compute the resampling ratio
            % required to achieve this in the pulse compressed time domain.
            if 0
              % ENABLE_FOR_DEBUG_OUTPUT_TIME_SAMPLING
              wf = 1;
              img = 1;
              rec = 1;
              wfs(wf).fs_raw = 250e6;
              wfs(wf).chirp_rate = 16e9/240e-6;
              wfs(wf).BW_window = [2.7e9 17.5e9];
              hdr.DDC_dec{img}(rec) = 5;
              
              wfs(wf).f0 = 2e9;
              wfs(wf).f1 = 18e9;
              wfs(wf).Tpd = 240e-6;
              BW = wfs(wf).f1-wfs(wf).f0;
              wfs(wf).chirp_rate = BW/wfs(wf).Tpd;
              hdr.t_ref{img}(rec) = 0e-6;
              tguard = 1e-6;
              td_max = 4e-6;
              td_min = 1e-6;
              tref = hdr.t_ref{img}(rec);
              hdr.t0_raw{img}(rec) = -tguard + min(td_max,tref);
              wfs(wf).Tadc_adjust = 0;
              fs_nyquist = max(wfs(wf).f0,wfs(wf).f1)*5;
              Mt_oversample = ceil(fs_nyquist/wfs(wf).fs_raw);
              fs_rf = Mt_oversample*wfs(wf).fs_raw;
              Nt_rf = round((wfs(wf).Tpd + max(td_max,tref) - min(td_min,tref) + 2*tguard)*fs_rf/Mt_oversample)*Mt_oversample;
              hdr.Nt{img}(rec) = Nt_rf / Mt_oversample;
            end
            Nt_desired = round(wfs(wf).fs_raw/wfs(wf).chirp_rate*diff(wfs(wf).BW_window)/2)*2;
            fs_raw_dec = wfs(wf).fs_raw ./ hdr.DDC_dec{img}(rec);
            Nt_raw_trim = round(fs_raw_dec/wfs(wf).chirp_rate*diff(wfs(wf).BW_window)/2)*2;
            [p,q] = rat(Nt_raw_trim*hdr.DDC_dec{img}(rec) / Nt_desired);
            % Create raw time domain window
            H_Nt = wfs(wf).ft_wind(Nt_raw_trim);
            % Create original raw time axis
            t0 = hdr.t0_raw{img}(rec) + wfs(wf).Tadc_adjust;
            dt_raw = 1/fs_raw_dec;
            time_raw_no_trim = (t0:dt_raw:t0+dt_raw*(hdr.Nt{img}(rec)-1)).';
            % Create RF frequency axis for minimum delay to surface expected
            f_rf = wfs(wf).f0 + wfs(wf).chirp_rate*(time_raw_no_trim - wfs(wf).td_mean);
            window_start_idx_norm = find(f_rf >= wfs(wf).BW_window(1),1);
            
            if 0
              % ENABLE_FOR_DEBUG_OUTPUT_TIME_SAMPLING
              df_before_resample = 1 / (Nt_raw_trim/wfs(wf).fs_raw*hdr.DDC_dec{img}(rec))
              df_after_resample = df_before_resample * p / q
              df_desired = 1 / (Nt_desired/wfs(wf).fs_raw)
              Nt_desired
              Nt_raw_trim
              p
              q
            end
            
            %% Determine the mapping of frequencies to time delays
            % =============================================================
            if 0
              % ENABLE_FOR_DEBUG_FREQ_MAP
              img = 1;
              rec = 1;
              hdr.Nt{img}(rec) = 10000;
              hdr.DDC_dec{img}(rec) = 3;
              wfs(wf).fs_raw = 100e6;
              hdr.nyquist_zone_signal{img}(rec) = 1;
              hdr.DDC_freq{img}(rec) = 95e6;
            end
            
            df_raw = wfs(wf).fs_raw/hdr.DDC_dec{img}(rec)/Nt_raw_trim;
            
            % nz: Nyquist zone containing signal spectrum (just renaming
            %   variable for convenience). The assumption is that the
            %   signal does not cross nyquist zones.
            nz = hdr.nyquist_zone_signal{img}(rec);
            
            % f_nz0: Lowest frequency in terms of ADC input frequency of the
            %   nyquist zone which contains the signal
            f_nz0 = wfs(wf).fs_raw * floor(nz/2);
            
            freq_raw =  f_nz0 + mod(hdr.DDC_freq{img}(rec) ...
              + df_raw*ifftshift(-floor(Nt_raw_trim/2):floor((Nt_raw_trim-1)/2)).', wfs(wf).fs_raw);
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
            if mod(length(unique_idxs),2)
              unique_idxs = unique_idxs(1:end-1);
            end
            
            freq_raw_unique = freq_raw_valid(unique_idxs);
            conjugate_unique = conjugate_bins(unique_idxs);
            % unique_idxs(~valid_bins);
            
            if 0
              % ENABLE_FOR_DEBUG_FREQ_MAP
              figure(1);
              clf;
              plot(freq_raw_valid,'.')
              hold on
              freq_raw_valid(conjugate_bins) = NaN;
              plot(freq_raw_valid,'r.')
              grid on
              xlabel('Bin (original order)');
              ylabel('Frequency (Hz)');
              legend('conj','-','location','best');
              
              figure(2);
              clf;
              plot(freq_raw_unique,'.')
              hold on
              tmp = freq_raw_unique;
              tmp(~conjugate_unique) = NaN;
              plot(tmp,'r.')
              grid on
              xlabel('Bin (reordered)');
              ylabel('Frequency (Hz)');
            end
          end
          
          % Check to see if reference time offset has changed since last record
          if freq_axes_changed ...
              || hdr.t_ref{img}(rec) ~= hdr.t_ref{img}(rec-1)
            
            % Convert IF frequency to time delay and account for reference
            % deramp time offset, hdr.t_ref
            time_rel = freq_raw_unique/wfs(wf).chirp_rate; % Relative time to t_ref
            time = time_rel + hdr.t_ref{img}(rec);
            
            % Ensure that start time is a multiple of dt
            dt = time(2)-time(1);
            time_correction = dt - mod(time(1),dt);
            time = time + time_correction;
            time_correction_IF = -time_correction*wfs(wf).chirp_rate;
            
            fc = sum(wfs(wf).BW_window)/2;
            Nt = length(time);
            T = Nt*dt;
            df = 1/T;
            freq = fc + df * ifftshift(-floor(Nt/2) : floor((Nt-1)/2)).';
            
            freq_axes_changed = false; % Reset state
          end
          
          % Get the start time for this record
          hdr.t0{img}(rec) = time(1);
          
          % Convert raw time into instantaneous frequency for the surface bin
          f_rf = wfs(wf).f0 + wfs(wf).chirp_rate*(time_raw_no_trim - hdr.surface(rec));
          if wfs(wf).BW_window(2) > max(f_rf)
            if ~BW_window_max_warning_printed
              BW_window_max_warning_printed = true;
              warning('BW_window (%g) is more than maximum measured RF frequency (%g) with surface twtt %g.', ...
                wfs(wf).BW_window(2), max(f_rf), hdr.surface(rec))
            end
            % Mark record as bad, but keep 2 bins to simplify later code
            hdr.Nt{img}(rec) = 2;
            data{img}(1:hdr.Nt{img}(rec),rec,wf_adc) = NaN;
            continue
          end
          if wfs(wf).BW_window(1) < min(f_rf)
            if ~BW_window_min_warning_printed
              BW_window_min_warning_printed = true;
              warning('BW_window (%g) is less than minimum measured RF frequency (%g) with surface twtt %g.', ...
                wfs(wf).BW_window(1), min(f_rf), hdr.surface(rec))
            end
            % Mark record as bad, but keep 2 bins to simplify later code
            hdr.Nt{img}(rec) = 2;
            data{img}(1:hdr.Nt{img}(rec),rec,wf_adc) = NaN;
            continue
          end
          
          % Create the window for the particular range line
          window_start_idx = find(f_rf >= wfs(wf).BW_window(1),1);
          H_idxs = window_start_idx : window_start_idx+Nt_raw_trim-1;
          if 0
            % ENABLE_FOR_DEBUG
            fprintf('window_start_idx: %d window_start_idx_norm: %d\n', ...
              window_start_idx, window_start_idx_norm);
          end
          
          % Window and Pulse compress
          tmp = data{img}(H_idxs,rec,wf_adc);
          if 0
            % Debug: Verify pulse compression window is correct
            clf;
            plot(lp(tmp));
            hold on;
            plot(lp(tmp .* H_Nt));
          end
          tmp = fft(tmp .* H_Nt);
          
          % Deskew of the residual video phase (not the standard because we
          % actually move the window to track the td)
          tmp_freq = (0:length(tmp)-1).';
          tmp = tmp .* exp(1i*2*pi*tmp_freq*(window_start_idx_norm-window_start_idx)/length(tmp));
          
          % Reorder result in case it is wrapped
          tmp = tmp(unique_idxs);
          % Some of the frequency bins are conjugated versions of the
          % signal
          tmp(conjugate_unique) = conj(tmp(conjugate_unique));
          
          % Complex baseband data (shifts by ~Tpd/2)
          tmp = ifftshift(fft(conj(tmp)));
          
          % Modulate the raw data to adjust the start time to always be a
          % multiple of wfs(wf).dt. Since we want this adjustment to be a
          % pure time shift and not introduce any phase shift in the other
          % domain, we make sure the phase is zero in the center of the
          % window: -time_raw(1+floor(Nt/2))
          tmp = tmp .* exp(1i*2*pi*freq*time_correction);
          
          % Return to time domain
          tmp = ifft(tmp);
          
          % Deskew of the residual video phase (second stage)
          tmp = tmp .* exp(-1i*pi*wfs(wf).chirp_rate*(time-wfs(wf).td_mean).^2);
          
          % Resample data so it aligns to constant time step
          if p~=q
            tmp = resample(tmp,p,q);
          end
          
          % Update the data matrix with the pulse compressed waveform
          hdr.Nt{img}(rec) = length(tmp);
          data{img}(1:hdr.Nt{img}(rec),rec,wf_adc) = tmp;
        end
        
        if 0
          % ENABLE_FOR_DEBUG
          figure(1); clf;
          Mt = 10;
          data_oversampled = interpft(data{img}(1:hdr.Nt{img}(rec),:,wf_adc), hdr.Nt{img}(rec)*Mt);
          [~,idx] = max(data_oversampled);
          time_oversampled = time(1) + dt/Mt* (0:length(time)*Mt-1).';
          plot((time_oversampled(idx).' - tds)/dt)
          grid on;
          xlabel('Record');
          ylabel('Time error (\Delta_t)');
          
          figure(2); clf;
          phase_sim = max(data{img}(1:hdr.Nt{img}(rec),:,wf_adc));
          plot(angle(phase_sim./phase_sim(1)),'+-');
          hold on
          fc_window = mean(wfs(wf).BW_window);
          phase_theory = exp(-1i*2*pi*fc_window*tds);
          plot(angle(phase_theory./phase_theory(1)),'.--');
          xlabel('Range bin');
          ylabel('Phase (rad)');
          legend('Simulated','Theory');
          grid on;
        end
        
        % Create a matrix of data with constant time rows, fill invalid samples with NaN
        if wf_adc == 1
          idx_start = min(round(hdr.t0{img}/dt));
          wfs(wf).Nt = max(round(hdr.t0{img}/dt) + hdr.Nt{img})-idx_start;
          hdr.time{img} = idx_start*dt + dt*(0:wfs(wf).Nt-1).';
          fc = sum(wfs(wf).BW_window)/2;
          T = wfs(wf).Nt*dt;
          df = 1/T;
          hdr.freq{img} = fc + df * ifftshift(-floor(wfs(wf).Nt/2) : floor((wfs(wf).Nt-1)/2)).';
        end
        % Method of copying to make this more efficient for very large
        % complex (real/imag) arrays. Lots of small matrix operations on
        % huge complex matrices is very slow in matlab. Real only matrices
        % are very fast though.
        blocks = round(linspace(1,size(data{img},2)+1,8)); blocks = unique(blocks);
        for block = 1:length(blocks)-1
          rlines = blocks(block) : blocks(block+1)-1;
          reD = real(data{img}(:,rlines,wf_adc));
          imD = imag(data{img}(:,rlines,wf_adc));
          for rec = 1:length(rlines)
            cur_idx_start = round(hdr.t0{img}(rlines(rec))/dt) - idx_start + 1;
            cur_idx_stop = round(hdr.t0{img}(rlines(rec))/dt) - idx_start + hdr.Nt{img}(rlines(rec));
            
            reD(cur_idx_start : cur_idx_stop,rec,wf_adc) = reD(1:hdr.Nt{img}(rlines(rec)),rec,wf_adc);
            reD(1:cur_idx_start-1,rec,wf_adc) = NaN;
            reD(cur_idx_stop+1 : wfs(wf).Nt,rec,wf_adc) = NaN;
          end
          for rec = 1:length(rlines)
            cur_idx_start = round(hdr.t0{img}(rlines(rec))/dt) - idx_start + 1;
            cur_idx_stop = round(hdr.t0{img}(rlines(rec))/dt) - idx_start + hdr.Nt{img}(rlines(rec));
            
            imD(cur_idx_start : cur_idx_stop,rec,wf_adc) = imD(1:hdr.Nt{img}(rlines(rec)),rec,wf_adc);
            imD(1:cur_idx_start-1,rec,wf_adc) = NaN;
            imD(cur_idx_stop+1 : wfs(wf).Nt,rec,wf_adc) = NaN;
          end
          data{img}(1:wfs(wf).Nt,rlines,wf_adc) = reD(1:wfs(wf).Nt,:) + 1i*imD(1:wfs(wf).Nt,:);
        end
        clear reD imD;
        
      elseif strcmpi(radar_type,'stepped')
        
      end
      
    else
      if wf_adc == 1
        if strcmpi(radar_type,'pulsed')
          hdr.time{img} = wfs(wf).time_raw;
          hdr.freq{img} = wfs(wf).freq_raw;
          
        elseif strcmpi(radar_type,'deramp')
          % Time axis is not valid if DDC or time offset changes
          hdr.time{img} = hdr.t0{img}(1) + 1/wfs(wf).fs_raw*(0:hdr.Nt{img}-1).';
          % Frequency is not valid
          df = wfs(wf).fs_raw / hdr.Nt{img};
          hdr.freq{img} = df*(0:hdr.Nt{img}-1).';
          
        elseif strcmpi(radar_type,'stepped')
        end
      end
    end
    
    %% Remove coherent noise for deramp radar
    % ===================================================================
    if strcmpi(radar_type,'deramp')
      
      if strcmpi(wfs(wf).coh_noise_method,'analysis')
        cdf_fn_dir = fileparts(ct_filename_out(param,wfs(wf).coh_noise_arg.fn, ''));
        cdf_fn = fullfile(cdf_fn_dir,sprintf('coh_noise_simp_%s_wf_%d_adc_%d.nc', param.day_seg, wf, adc));
        
        finfo = ncinfo(cdf_fn);

        noise = [];
        coh_noise_date_str = ncread(cdf_fn,'datestr');
        hdr.custom.coh_noise(1:length(coh_noise_date_str),1,img,wf_adc) = coh_noise_date_str;
        noise.start_bin = ncread(cdf_fn,'start_bin');
        noise.dt = ncread(cdf_fn,'dt');
        noise.dft_freqs = ncread(cdf_fn,'dft_freqs');
        noise.recs = ncread(cdf_fn,'recs');
        noise.gps_time = ncread(cdf_fn,'gps_time');
        noise.Nx = length(noise.gps_time);
        
        tmp = netcdf_to_mat(cdf_fn,[],'^param_analysis');
        noise.param_analysis = tmp.param_analysis;
        
        dt = hdr.time{img}(2)-hdr.time{img}(1);
        if abs(noise.dt-dt)/dt > 1e-6
          error('There is a fast-time sample interval discrepancy between the current processing settings (%g) and those used to generate the coherent noise file (%g).', dt, noise.dt);
        end
        % Adjust coherent noise start_bin for changes in t_ref relative to
        % when the coherent noise was loaded and estimated.
        delta_t_ref_bin = (noise.param_analysis.radar.wfs(wf).t_ref - wfs(wf).t_ref)/dt;
        if abs(round(delta_t_ref_bin)-delta_t_ref_bin) > 1e-3
          error('There is a fast-time reference time delay discrepancy between the current processing settings (%g) and those used to generate the coherent noise file (%g). Changes in t_ref require recreating the coherent noise file or the changes must be a multiple of a time bin (%g).', ...
            wfs(wf).t_ref, noise.param_analysis.radar.wfs(wf).t_ref, dt);
        end
        start_bin = round(hdr.time{img}(1)/dt) + delta_t_ref_bin;
        
        % Nt by Nx_dft matrix (we grab a subset of the Nt samples)
        noise.dft = ncread(cdf_fn,'dftI',[start_bin-noise.start_bin+1 1],[wfs(wf).Nt inf]) ...
          + 1i*ncread(cdf_fn,'dftQ',[start_bin-noise.start_bin+1 1],[wfs(wf).Nt inf]);
        
        % Adjust coherent noise dft for changes in adc_gains relative to
        % when the coherent noise was loaded and estimated.
        noise.dft = noise.dft * wfs(wf).adc_gains(adc) ./ noise.param_analysis.radar.wfs(wf).adc_gains(adc);
        
        recs = interp1(noise.gps_time, noise.recs, hdr.gps_time);
        
        if 0
          imagesc(lp( bsxfun(@minus, data{1}(1:wfs(wf).Nt,:), mean(data{1}(1:wfs(wf).Nt,:),2) )  ))
          figure(1); clf;
          plot(lp(mean(data{1}(1:wfs(wf).Nt,:),2)))
          hold on
          plot(lp(noise.dft))
        end
        
        for dft_idx = 1:length(noise.dft_freqs)
          % mf: matched filter
          % noise.dft(bin,dft_idx): Coefficient for the matched filter
          mf = exp(1i*2*pi/noise.Nx*noise.dft_freqs(dft_idx) .* recs);
          for bin = 1:wfs(wf).Nt
            data{img}(bin,:,wf_adc) = data{img}(bin,:,wf_adc) ...
              - noise.dft(bin,dft_idx) * mf;
          end
        end
        
      elseif strcmpi(wfs(wf).coh_noise_method,'estimated')
        % Apply coherent noise methods that require estimates derived now
        
        if wfs(wf).coh_noise_arg.DC_remove_en
          data{img}(1:wfs(wf).Nt,:,wf_adc) = bsxfun(@minus, data{img}(1:wfs(wf).Nt,:,wf_adc), ...
            mean(data{img}(1:wfs(wf).Nt,:,wf_adc),2));
        end
        
        if length(wfs(wf).coh_noise_arg.B_coh_noise) > 1
          if length(wfs(wf).coh_noise_arg.A_coh_noise) > 1
            % Use filtfilt (feedback)
            data{img}(1:wfs(wf).Nt,:,wf_adc) = single(filtfilt(wfs(wf).coh_noise_arg.B_coh_noise, ...
              wfs(wf).coh_noise_arg.A_coh_noise, double(data{img}(1:wfs(wf).Nt,:,wf_adc).'))).';
          else
            % Use fir_dec (no feedback)
            data{img}(1:wfs(wf).Nt,:,wf_adc) = fir_dec(data{img}(1:wfs(wf).Nt,:,wf_adc),wfs(wf).coh_noise_arg.B_coh_noise,1);
          end
        end
        
      end
    end
    
    %% Deconvolution
    % ===================================================================
    
    if wfs(wf).deconv.en
      deconv_fn = fullfile(fileparts(ct_filename_out(param,wfs(wf).deconv.fn, '')), ...
        sprintf('deconv_%s_wf_%d_adc_%d.mat',param.day_seg, wf, adc));
      deconv = load(deconv_fn);
      deconv_date_str = deconv.param_collate_deconv_final.sw_version.cur_date_time;
      hdr.custom.deconv(1:length(deconv_date_str),1,img,wf_adc) = deconv_date_str;
      
      deconv_map_idxs = interp1(deconv.map_gps_time,deconv.map_idxs,hdr.gps_time,'nearest','extrap');
      max_score = interp1(deconv.map_gps_time,deconv.max_score,hdr.gps_time,'nearest','extrap');
      
      unique_idxs = unique(deconv_map_idxs);
      
      if wf_adc == 1
        % Prepare variables that are constant for each deconvolution
        % waveform that will be applied
        cmd = deconv.param_collate_deconv.analysis.cmd{deconv.param_collate_deconv.collate_deconv.cmd_idx};
        freq = fftshift(hdr.freq{img});
      fc = hdr.freq{img}(1);
        
        % Prepare variables to baseband data to new center frequency (in
        % case the deconvolution filter subbands)
        deconv_fc = (cmd.f0+cmd.f1)/2;
        df = hdr.freq{img}(2)-hdr.freq{img}(1);
        BW = df * wfs(wf).Nt;
        deconv_dfc = deconv_fc - fc;
        hdr.freq{img} = mod(hdr.freq{img} + deconv_dfc-wfs(wf).BW_window(1), BW)+wfs(wf).BW_window(1);
        
      elseif abs(deconv_fc - (cmd.f0+cmd.f1)/2)/deconv_fc > 1e-6
        error('Deconvolution center frequency must be the same for all wf-adc pairs in the image. Was %g and is now %g.', deconv_fc, (cmd.f0+cmd.f1)/2);
      end
      
      for unique_idxs_idx = 1:length(unique_idxs)
        % deconv_mask: Create logical mask corresponding to range lines that use this deconv waveform
        deconv_map_idx = unique_idxs(unique_idxs_idx);
        deconv_mask = deconv_map_idx == deconv_map_idxs;
        % deconv_mask = deconv_map_idx == deconv_map_idxs ...
        %   & max_score > deconv.param_collate_deconv_final.collate_deconv.min_score;
        
        if wfs(wf).Nt <= 2 || ~any(deconv_mask)
          % Range lines are bad (Nt <= 2), or no matching range lines that
          % have a good enough score to justify deconvolution.
          continue;
        end
        
        % Get the reference function
        h_nonnegative = deconv.ref_nonnegative{deconv_map_idx};
        h_negative = deconv.ref_negative{deconv_map_idx};
        h_mult_factor = deconv.ref_mult_factor(deconv_map_idx);
        
        % Adjust deconvolution signal to match sample rline
        h_filled = [h_nonnegative; zeros(wfs(wf).Nt-length(h_nonnegative)-length(h_negative),1); h_negative];
        
        % Is dt different? Error
        dt = hdr.time{img}(2)-hdr.time{img}(1);
        if abs(deconv.dt-dt)/dt > 1e-6
          error('There is fast-time sample interval discrepancy between the current processing settings (%g) and those used to generate the deconvolution file (%g).', dt, deconv.dt);
        end
        
        % Is fc different? Multiply time domain by exp(1i*2*pi*dfc*deconv_time)
        dfc = fc - deconv.fc(deconv_map_idx);
        if dfc/fc > 1e-6
          deconv_time = t0 + dt*(0:Nt-1).';
          h_filled = h_filled .* exp(1i*2*pi*dfc*deconv_time);
        end
        deconv_LO = exp(-1i*2*pi*(dfc+deconv_dfc) * hdr.time{img});
        
        % Take FFT of deconvolution impulse response
        h_filled = fft(h_filled);
        
        % Create inverse filter relative to window
        Nt_shorten = find(cmd.f0 <= freq,1);
        Nt_shorten(2) = length(freq) - find(cmd.f1 >= freq,1,'last');
        Nt_Hwind = wfs(wf).Nt - sum(Nt_shorten);
        Hwind = deconv.ref_window(Nt_Hwind);
        Hwind_filled = ifftshift([zeros(Nt_shorten(1),1); Hwind; zeros(Nt_shorten(end),1)]);
        h_filled_inverse = Hwind_filled ./ h_filled;
        
        % Normalize deconvolution
        h_filled_inverse = h_filled_inverse * h_mult_factor;
        
        % Apply deconvolution filter
        deconv_mask_idxs = find(deconv_mask);
        blocks = round(linspace(1,length(deconv_mask_idxs)+1,8)); blocks = unique(blocks);
        for block = 1:length(blocks)-1
          rlines = deconv_mask_idxs(blocks(block) : blocks(block+1)-1);
          % Matched filter
          data{img}(1:wfs(wf).Nt,rlines,wf_adc) = ifft(bsxfun(@times, fft(data{img}(1:wfs(wf).Nt,rlines,wf_adc)), h_filled_inverse));
          % Down conversion to new deconvolution center frequency
          data{img}(1:wfs(wf).Nt,rlines,wf_adc) = bsxfun(@times, data{img}(1:wfs(wf).Nt,rlines,wf_adc), deconv_LO);
        end
        
      end
      
    end
    
  end
  
  if param.load.pulse_comp == 1
    data{img} = data{img}(1:wfs(wf).Nt,:,:);
  end
end
