function [success] = analysis_combine_task(param)
% [success] = analysis_combine_task(param)
%
% Cluster task for combining results from analysis.
%
% param: struct controlling the processing
%
% success = boolean which is true when the function executes properly
%   if a task fails before it can return success, then success will be
%   empty
%
% Author: John Paden
%
% See also analysis.m

% Load records file
records_fn = ct_filename_support(param,'','records');
records = load(records_fn);
% Apply presumming
if param.analysis.presums > 1
  records.lat = fir_dec(records.lat,param.analysis.presums);
  records.lon = fir_dec(records.lon,param.analysis.presums);
  records.elev = fir_dec(records.elev,param.analysis.presums);
  records.roll = fir_dec(records.roll,param.analysis.presums);
  records.pitch = fir_dec(records.pitch,param.analysis.presums);
  records.heading = fir_dec(records.heading,param.analysis.presums);
  records.gps_time = fir_dec(records.gps_time,param.analysis.presums);
  records.surface = fir_dec(records.surface,param.analysis.presums);
end

[~,~,radar_name] = ct_output_dir(param.radar_name);

% Break records in segment into blocks
breaks = 1:param.analysis.block_size:length(records.gps_time);

% If the last block is less than half the desired block size, then combine
% with earlier block if possible
if length(records.gps_time)-breaks(end) < param.analysis.block_size/2 ...
    && length(breaks) > 1
  breaks = breaks(1:end-1);
end

%% Loop through all given commands from 'analysis'
for cmd_idx = 1:length(param.analysis.cmd)
  cmd = param.analysis.cmd{cmd_idx};
  if ~cmd.en
    continue;
  end
  
  if strcmpi(cmd.method,{'saturation'})
    %% Saturation check
    % ===================================================================
    % ===================================================================
    
    for img = 1:length(param.analysis.imgs)
      max_rline = [];
      max_waveform = [];
      gps_time = [];
      max_val_gps_time = [];
      max_val_gps_time_adc = [];
      
      for break_idx = 1:length(breaks)
        rec_load_start = breaks(break_idx);
        
        if break_idx == length(breaks)
          rec_load_stop = length(records.gps_time);
        else
          rec_load_stop = rec_load_start+param.analysis.block_size-1;
        end
        
        % =====================================================================
        % Prepare task inputs
        % =====================================================================
        cur_recs = [rec_load_start rec_load_stop];
        actual_cur_recs = [(cur_recs(1)-1)*param.analysis.presums+1, ...
          cur_recs(end)*param.analysis.presums];
        
        out_fn = fullfile(ct_filename_out(param, param.analysis.out_path), ...
          sprintf('saturation_img_%02d_%d_%d.mat',img,actual_cur_recs));
        
        satur = load(out_fn);
        
        max_rline = cat(2,max_rline,satur.max_rline);
        max_waveform = cat(2,max_waveform,satur.max_waveform);
        gps_time = cat(2,gps_time,satur.gps_time);
        max_val_gps_time = cat(2,max_val_gps_time,satur.max_val_gps_time);
        if isinteger(satur.max_val_gps_time)
          max_val_gps_time_adc = cat(2,max_val_gps_time_adc,satur.max_val_gps_time_adc);
        end
        
        satur.max_rline = max_rline;
        satur.max_waveform = max_waveform;
        satur.gps_time = gps_time;
        satur.max_val_gps_time = max_val_gps_time;
        satur.max_val_gps_time_adc = max_val_gps_time_adc;
        
        out_fn_dir = fileparts(out_fn);
        out_segment_fn_dir = fileparts(out_fn_dir);
        out_segment_fn = fullfile(out_segment_fn_dir,sprintf('saturation_%s_img_%02d.mat', param.day_seg, img));
        fprintf('Saving output %s (%s)\n', out_segment_fn, datestr(now));
        save(out_segment_fn,'-v7.3','-struct','satur');
      end
    end
    
    
  elseif strcmpi(cmd.method,{'specular'})
    %% Specular Analysis for Deconvolution
    % ===================================================================
    % ===================================================================
    for img = 1:length(param.analysis.imgs)
      for wf_adc = 1:size(param.analysis.imgs{img},1)
        wf = param.analysis.imgs{1}(wf_adc,1);
        adc = param.analysis.imgs{1}(wf_adc,2);
        
        gps_time = [];
        lat = [];
        lon = [];
        elev = [];
        roll = [];
        pitch = [];
        heading = [];
        peakiness = [];
        deconv_gps_time = [];
        deconv_mean = {};
        deconv_std = {};
        deconv_sample = {};
        deconv_freq = {};
        deconv_t0 = [];
        deconv_twtt = [];
        deconv_forced = [];
        deconv_DDC_Mt = [];
        for break_idx = 1:length(breaks)
          rec_load_start = breaks(break_idx);
          
          if break_idx == length(breaks)
            rec_load_stop = length(records.gps_time);
          else
            rec_load_stop = rec_load_start+param.analysis.block_size-1;
          end
          
          %% Specular: Load task output and concatenate
          % ===============================================================
          cur_recs = [rec_load_start rec_load_stop];
          actual_cur_recs = [(cur_recs(1)-1)*param.analysis.presums+1, ...
            cur_recs(end)*param.analysis.presums];
          
          out_fn = fullfile(ct_filename_out(param, param.analysis.out_path), ...
            sprintf('specular_wf_%d_adc_%d_%d_%d.mat',wf,adc,actual_cur_recs));
          
          spec = load(out_fn);
          
          fprintf('%s\n',datestr(now))
          
          % Concatenate
          gps_time(end+(1:numel(spec.gps_time))) = spec.gps_time;
          lat(end+(1:numel(spec.lat))) = spec.lat;
          lon(end+(1:numel(spec.lon))) = spec.lon;
          elev(end+(1:numel(spec.elev))) = spec.elev;
          roll(end+(1:numel(spec.roll))) = spec.roll;
          pitch(end+(1:numel(spec.pitch))) = spec.pitch;
          heading(end+(1:numel(spec.heading))) = spec.heading;
          peakiness(end+(1:numel(spec.peakiness))) = spec.peakiness;
          
          deconv_gps_time(end+(1:numel(spec.deconv_gps_time))) = spec.deconv_gps_time;
          deconv_twtt(end+(1:numel(spec.deconv_twtt))) = spec.deconv_twtt;
          deconv_DDC_Mt(end+(1:numel(spec.deconv_DDC_Mt))) = spec.deconv_DDC_Mt;
          deconv_forced(end+(1:numel(spec.deconv_forced))) = spec.deconv_forced;
          
          deconv_mean(end+(1:numel(spec.deconv_mean))) = spec.deconv_mean;
          deconv_std(end+(1:numel(spec.deconv_std))) = spec.deconv_std;
          deconv_sample(end+(1:numel(spec.deconv_sample))) = spec.deconv_sample;
          for idx = 1:numel(spec.deconv_gps_time)
            deconv_freq{end+1} = spec.freq;
            deconv_t0(end+1) = spec.time(1);
          end
        end
        
        %% Specular: Store concatenated output
        % =================================================================
        spec.gps_time = gps_time;
        spec.lat = lat;
        spec.lon = lon;
        spec.elev = elev;
        spec.roll = roll;
        spec.pitch = pitch;
        spec.heading = heading;
        spec.peakiness = peakiness;
        spec.deconv_gps_time = deconv_gps_time;
        spec.deconv_mean = deconv_mean;
        spec.deconv_std = deconv_std;
        spec.deconv_sample = deconv_sample;
        spec.deconv_freq = deconv_freq;
        spec.deconv_t0 = deconv_t0;
        spec.deconv_twtt = deconv_twtt;
        spec.deconv_forced = deconv_forced;
        spec.deconv_DDC_Mt = deconv_DDC_Mt;
        spec = rmfield(spec, {'time','freq'});
        out_fn_dir = fileparts(out_fn);
        out_segment_fn_dir = fileparts(out_fn_dir);
        out_segment_fn = fullfile(out_segment_fn_dir,sprintf('specular_%s_wf_%d_adc_%d.mat', param.day_seg, wf, adc));
        fprintf('Saving output %s (%s)\n', out_segment_fn, datestr(now));
        save(out_segment_fn,'-v7.3','-struct','spec');
      end
    end
    
    
  elseif strcmpi(cmd.method,{'coh_noise'})
    %% Coherent Noise Analysis
    % ===================================================================
    % ===================================================================
    for img = 1:length(param.analysis.imgs)
      
      for wf_adc_idx = cmd.wf_adc_idxs{img}(:).'
        wf = param.analysis.imgs{1}(wf_adc_idx,1);
        adc = param.analysis.imgs{1}(wf_adc_idx,2);
        
        %% Loop through all the coherent noise tracker files and combine
        % =====================================================================
        gps_time = [];
        lat = [];
        lon = [];
        elev = [];
        roll = [];
        pitch = [];
        heading = [];
        nyquist_zone = [];
        coh_ave = [];
        coh_ave_samples = [];
        doppler_concat = [];
        for break_idx = 1:length(breaks)
          rec_load_start = breaks(break_idx);
          
          if break_idx == length(breaks)
            rec_load_stop = length(records.gps_time);
          else
            rec_load_stop = rec_load_start+param.analysis.block_size-1;
          end
          
          % =====================================================================
          % Prepare task inputs
          % =====================================================================
          cur_recs = [rec_load_start rec_load_stop];
          actual_cur_recs = [(cur_recs(1)-1)*param.analysis.presums+1, ...
            cur_recs(end)*param.analysis.presums];
          
          out_fn = fullfile(ct_filename_out(param, param.analysis.out_path), ...
            sprintf('coh_noise_wf_%d_adc_%d_%d_%d.mat',wf,adc,actual_cur_recs));
          
          noise = load(out_fn);
          if strcmpi(radar_name,'fmcw')
            if size(noise.coh_ave,1) ~= num_samples
              warning('A BAD RESAMPLING METHOD IS BEING APPLIED TO THE DATA AND LIKELY TO PRODUCE POOR RESULTS. THIS CODE NEEDS TO BE REPLACED');
              if any(any(isnan(noise.coh_ave))) || any(any(isnan(noise.coh_ave_samples)))
                warning('NaN found in noise.coh_ave or noise.coh_ave_samples')
              end
              noise.coh_ave = interp1([1:size(noise.coh_ave,1)],noise.coh_ave,linspace(1,size(noise.coh_ave,1),num_samples));
              noise.coh_ave_samples = interp1([1:size(noise.coh_ave_samples,1)],noise.coh_ave_samples,linspace(1,size(noise.coh_ave_samples,1),num_samples));
            end
          end
          
          gps_time = cat(2,gps_time,noise.gps_time);
          lat = cat(2,lat,noise.lat);
          lon = cat(2,lon,noise.lon);
          elev = cat(2,elev,noise.elev);
          roll = cat(2,roll,noise.roll);
          pitch = cat(2,pitch,noise.pitch);
          heading = cat(2,heading,noise.heading);
          nyquist_zone = cat(2,nyquist_zone,noise.nyquist_zone);
          coh_ave = cat(2,coh_ave,noise.coh_ave);
          coh_ave_samples = cat(2,coh_ave_samples,noise.coh_ave_samples);
          noise.doppler = reshape(noise.doppler,[numel(noise.doppler) 1]);
          if break_idx > 1 && size(noise.doppler,1) ~= size(doppler_concat,1)
            % Block was a different size than other Doppler spectrums, re-sample
            % so that it can be stored in the output matrix
            noise.doppler = interp1(0:numel(noise.doppler)-1,noise.doppler,linspace(0,numel(noise.doppler)-1,size(doppler_concat,1)).');
          end
          doppler_concat = cat(2,doppler_concat,noise.doppler);
          
        end
        
        noise.gps_time = gps_time;
        noise.lat = lat;
        noise.lon = lon;
        noise.elev = elev;
        noise.roll = roll;
        noise.pitch = pitch;
        noise.heading = heading;
        noise.nyquist_zone = nyquist_zone;
        noise.coh_ave = coh_ave;
        noise.coh_ave_samples = coh_ave_samples;
        noise.doppler = doppler_concat;
        out_fn_dir = fileparts(out_fn);
        out_segment_fn_dir = fileparts(out_fn_dir);
        out_segment_fn = fullfile(out_segment_fn_dir,sprintf('coh_noise_%s_wf_%d_adc_%d.mat', param.day_seg, wf, adc));
        fprintf('Saving output %s (%s)\n', out_segment_fn, datestr(now));
        save(out_segment_fn,'-v7.3','-struct','noise'); % Use HDF because of the large file size
      end
    end
    
    
  elseif strcmpi(cmd.method,{'waveform'})
    %% Waveform extraction
    % ===================================================================
    % ===================================================================
    
    %% Loop through all the surface tracker files and combine
    % =====================================================================
    for img = 1:length(param.analysis.imgs)
      gps_time = [];
      lat = [];
      lon = [];
      elev = [];
      roll = [];
      pitch = [];
      heading = [];
      surf_vals = [];
      surf_bins = [];
      for break_idx = 1:length(breaks)
        rec_load_start = breaks(break_idx);
        
        if break_idx == length(breaks)
          rec_load_stop = length(records.gps_time);
        else
          rec_load_stop = rec_load_start+param.analysis.block_size-1;
        end
        
        % =====================================================================
        % Prepare task inputs
        % =====================================================================
        cur_recs = [rec_load_start rec_load_stop];
        actual_cur_recs = [(cur_recs(1)-1)*param.analysis.presums+1, ...
          cur_recs(end)*param.analysis.presums];    
        
        out_fn = fullfile(ct_filename_out(param, param.analysis.out_path), ...
          sprintf('surf_img_%02d_%d_%d.mat',img,actual_cur_recs));
        
        surf = load(out_fn);
        
        gps_time = cat(2,gps_time,surf.gps_time);
        lat = cat(2,lat,surf.lat);
        lon = cat(2,lon,surf.lon);
        elev = cat(2,elev,surf.elev);
        roll = cat(2,roll,surf.roll);
        pitch = cat(2,pitch,surf.pitch);
        heading = cat(2,heading,surf.heading);
        surf_vals = cat(2,surf_vals,surf.surf_vals);
        surf_bins = cat(2,surf_bins,surf.surf_bins);
      end
      
      surf.gps_time = gps_time;
      surf.lat = lat;
      surf.lon = lon;
      surf.elev = elev;
      surf.roll = roll;
      surf.pitch = pitch;
      surf.heading = heading;
      surf.surf_vals = surf_vals;
      surf.surf_bins = surf_bins;
      
      out_fn_dir = fileparts(out_fn);
      out_segment_fn_dir = fileparts(out_fn_dir);
      out_segment_fn = fullfile(out_segment_fn_dir,sprintf('surf_%s_img_%02d.mat', param.day_seg, img));
      fprintf('Saving output %s (%s)\n', out_segment_fn, datestr(now));
      save(out_segment_fn,'-v7.3','-struct','surf');
    end
    
    
  elseif strcmpi(cmd.method,{'statistics'})
    %% Statistical analysis
    % ===================================================================
    % ===================================================================
    
    %% Loop through all the power files and combine
    % =====================================================================
    for img = 1:length(param.analysis.imgs)
      gps_time = [];
      lat = [];
      lon = [];
      elev = [];
      roll = [];
      pitch = [];
      heading = [];
      power_vals = [];
      power_bins = [];
      for break_idx = 1:length(breaks)
        rec_load_start = breaks(break_idx);
        
        if break_idx == length(breaks)
          rec_load_stop = length(records.gps_time);
        else
          rec_load_stop = rec_load_start+param.analysis.block_size-1;
        end
        
        % =====================================================================
        % Prepare task inputs
        % =====================================================================
        cur_recs = [rec_load_start rec_load_stop];
        actual_cur_recs = [(cur_recs(1)-1)*param.analysis.presums+1, ...
          cur_recs(end)*param.analysis.presums];
        
        out_fn = fullfile(ct_filename_out(param, param.analysis.out_path), ...
          sprintf('statistics_img_%02d_%d_%d.mat',img,actual_cur_recs));
        
        power = load(out_fn);
        
        gps_time = cat(2,gps_time,power.gps_time);
        lat = cat(2,lat,power.lat);
        lon = cat(2,lon,power.lon);
        elev = cat(2,elev,power.elev);
        roll = cat(2,roll,power.roll);
        pitch = cat(2,pitch,power.pitch);
        heading = cat(2,heading,power.heading);
        power_vals = cat(2,power_vals,power.power_vals);
        power_bins = cat(2,power_bins,power.power_bins);
      end
      
      power.gps_time = gps_time;
      power.lat = lat;
      power.lon = lon;
      power.elev = elev;
      power.roll = roll;
      power.pitch = pitch;
      power.heading = heading;
      power.power_vals = power_vals;
      power.power_bins = power_bins;
      
      out_fn_dir = fileparts(out_fn);
      out_segment_fn_dir = fileparts(out_fn_dir);
      out_segment_fn = fullfile(out_segment_fn_dir,sprintf('power_%s_img_%02d.mat', param.day_seg, img));
      fprintf('Saving output %s (%s)\n', out_segment_fn, datestr(now));
      save(out_segment_fn,'-v7.3','-struct','power');
    end
    
    %% Loop through all the psd (power spectral density) files and combine
    % =====================================================================
    for img = 1:length(param.analysis.imgs)
      gps_time = [];
      lat = [];
      lon = [];
      elev = [];
      roll = [];
      pitch = [];
      heading = [];
      psd_vals = [];
      psd_bins = [];
      psd_mean = [];
      psd_Rnn = [];
      for break_idx = 1:length(breaks)
        rec_load_start = breaks(break_idx);
        
        if break_idx == length(breaks)
          rec_load_stop = length(records.gps_time);
        else
          rec_load_stop = rec_load_start+param.analysis.block_size-1;
        end
        
        % =====================================================================
        % Prepare task inputs
        % =====================================================================
        cur_recs = [rec_load_start rec_load_stop];
        
        out_fn = fullfile(ct_filename_out(param, ...
          param.analysis.out_path, 'CSARP_noise'), ...
          sprintf('psd_img_%02d_%d_%d.mat', img, cur_recs(1),cur_recs(end)));
        
        psd = load(out_fn);
        
        gps_time = cat(2,gps_time,psd.gps_time);
        lat = cat(2,lat,psd.lat);
        lon = cat(2,lon,psd.lon);
        elev = cat(2,elev,psd.elev);
        roll = cat(2,roll,psd.roll);
        pitch = cat(2,pitch,psd.pitch);
        heading = cat(2,heading,psd.heading);
        psd_vals = cat(2,psd_vals,psd.psd_vals);
        psd_bins = cat(2,psd_bins,psd.psd_bins);
        psd_mean = cat(2,psd_mean,psd.psd_mean);
        psd_Rnn = cat(2,psd_Rnn,psd.psd_Rnn);
      end
      
      psd.gps_time = gps_time;
      psd.lat = lat;
      psd.lon = lon;
      psd.elev = elev;
      psd.roll = roll;
      psd.pitch = pitch;
      psd.heading = heading;
      psd.psd_vals = psd_vals;
      psd.psd_bins = psd_bins;
      psd.psd_mean = psd_mean;
      psd.psd_Rnn = psd_Rnn;
      
      out_fn_dir = fileparts(out_fn);
      out_segment_fn_dir = fileparts(out_fn_dir);
      out_segment_fn = fullfile(out_segment_fn_dir,sprintf('psd_%s_img_%02d.mat', param.day_seg, img));
      fprintf('Saving output %s (%s)\n', out_segment_fn, datestr(now));
      save(out_segment_fn,'-v7.3','-struct','psd');
    end
    
  end
end

fprintf('%s done %s\n', mfilename, datestr(now));

success = true;

return;
