function [param] = load_icards_data(param,whole_param)

if ~isfield(param.proc,'raw_data')
  param.proc.raw_data = false;
end
if ~isfield(param.load,'wf_adc_comb')
  param.load.wf_adc_comb.en = 0;
end
sample_size=2;      %the file type of icards sample is "int16"---qishi
global g_data;
wfs = param.wfs;

physical_constants;

accum = build_img_load_struct(param.load.imgs, param.load.adcs);

% ===================================================================
% Preallocate data matrix
%====================================================================
total_rec = length(param.load.file_idx{1});
Nx = floor(total_rec/param.proc.presums);
if ~iscell(g_data)
  g_data = cell(size(param.load.imgs));
end
for img_idx = 1:length(param.load.imgs)
  wf = param.load.imgs{img_idx}(1,1);
  if param.load.wf_adc_comb.en
    Nt = param.load.wf_adc_comb.Nt;
  else
    Nt = param.wfs(wf).Nt;
  end
  if param.proc.combine_rx
    Nc = 1;
  else
    Nc = size(param.load.imgs{img_idx},1);
  end
  if any([size(g_data{img_idx},1) size(g_data{img_idx},2) size(g_data{img_idx},3)] ~= [Nt Nx Nc])
    g_data{img_idx} = zeros(Nt,Nx,Nc,'single');
  elseif param.proc.combine_rx
    g_data{img_idx}(:) = 0;
  end
end

% ===================================================================
% Load data
% ===================================================================
adc = 1;%only 1 for icards
adc_idx=1;
num_accum = 0;
out_idx = 0;
rec = 1;

if strcmpi(whole_param.radar.icards.data_type,'incoherent') %read coherent data
  if size(icards_get_data(param.load.filenames{1}{end},2),1)~=0%some coherent data in 1997's segments
    whole_param.radar.icards.data_type='coherent';
    warning('%s has coherent data, pay attention !\n',whole_param.day_seg);
  end
end
 
if strcmpi(whole_param.radar.icards.data_type,'coherent') %read coherent data
    while rec < total_rec;
      % Get the filename
      fn_idx = param.load.file_idx{adc_idx}(rec);
      fn = param.load.filenames{adc_idx}{fn_idx};

      % Get number of records to load
      num_rec = find(param.load.file_idx{adc_idx}(rec:total_rec ) > fn_idx, 1) - 1;

      % Check to see if we crossed a file boundary
      if isempty(num_rec)
        % Not crossing a file boundary
        num_rec = length(rec:total_rec);
        boundary_crossing = false;
      else
        % Crossing a file boundary
        boundary_crossing = true;
      end

      % ===============================================================
      % Load the records
      % ===============================================================

      fprintf('Loading data from %s\n', fn);
      [fid,msg] = fopen(fn, 'r');
      if fid<0
        error('File open failed (%s)\n%s',fn, msg);
      end
      file_record_length=size(icards_get_data(fn,2),2);

      for rec = rec:rec+num_rec-1

        % Load data record from file
        fseek(fid, param.load.offset{1}(rec), 'bof');
        rec_data_I = [fread(fid, param.load.rec_data_size/sample_size, 'int16')];%read I channel
        fseek(fid, param.load.offset{1}(rec)+param.load.rec_data_size*file_record_length+12, 'bof');
        rec_data_Q = [fread(fid, param.load.rec_data_size/sample_size, 'int16')];%read Q Channel
        if whole_param.radar.icards.IQ_flip
          rec_data = double(rec_data_Q+1i*rec_data_I);%a full sample I+jQ
        else
          rec_data = double(rec_data_I+1i*rec_data_Q);%a full sample I+jQ
        end

        % ===============================================================
        % Process record
        % ===============================================================
        accum_idx = 1; %icards just has one waveform
        if num_accum == 0
          accum.data{accum_idx} = rec_data;
        else
          accum.data{accum_idx} = accum.data{accum_idx} + rec_data;
        end

        num_accum = num_accum + 1;
        if num_accum >= param.proc.presums
          num_accum = 0;
          out_idx = out_idx + 1;

          wf = accum(adc).wf(accum_idx);
          img_idx = accum(adc).img_idx(accum_idx);
          wf_adc_idx = accum(adc).wf_adc_idx(accum_idx);

          % Apply channel compensation
          chan_equal = 10.^(param.radar.wfs(wf).chan_equal_dB(param.radar.wfs(wf).rx_paths(adc))/20) ...
            .* exp(j*param.radar.wfs(wf).chan_equal_deg(param.radar.wfs(wf).rx_paths(adc))/180*pi);
          accum(adc).data{accum_idx} = accum(adc).data{accum_idx}/chan_equal;
          accum(adc).data{accum_idx} = accum(adc).data{accum_idx}/wfs(wf).adc_gains(adc);
          if param.proc.pulse_comp
            % ===========================================================
            % Do pulse compression
            % Apply blank (only should enable if sidelobe problems present)
            if ~isempty(wfs(wf).blank)
              % accum(adc).data{accum_idx}(wfs(wf).time_raw>wfs(wf).blank(1) & wfs(wf).time_raw<wfs(wf).blank(2)) = 0;
              accum(adc).data{accum_idx}(wfs(wf).time_raw-param.radar.wfs(wf).Tsys(adc) <= param.surface(rec) + wfs(wf).blank) = 0;
            end
            % Apply matched filter
            % Zero pad front: (the standard)
            accum(adc).data{accum_idx} = fft([zeros(wfs(wf).pad_length,1); accum(adc).data{accum_idx}], wfs(wf).Nt_pc);
            % Zero pad end: (debug only)
            %accum(adc).data{accum_idx} = fft(accum(adc).data{accum_idx}, wfs(wf).Nt_pc);
            
            % Apply matched filter and transform back to time domain
            accum(adc).data{accum_idx} = ifft(accum(adc).data{accum_idx} .* wfs(wf).ref{adc});
            
            if param.proc.ft_dec
              % Digital down conversion and decimation
              accum(adc).data{accum_idx} = accum(adc).data{accum_idx}.*exp(-1i*2*pi*wfs(wf).fc*wfs(wf).time_raw);
              accum(adc).data{accum_idx} = resample(double(accum(adc).data{accum_idx}), param.wfs(1).ft_dec(1), param.wfs(1).ft_dec(2));
            end
            
          elseif param.proc.ft_dec
            accum(adc).data{accum_idx} = fft(accum(adc).data{accum_idx},wfs(wf).Nt_raw);
            accum(adc).data{accum_idx} = ifft(accum(adc).data{accum_idx}(wfs(wf).freq_inds));
            accum(adc).data{accum_idx} = accum(adc).data{accum_idx}.*exp(-1i*2*pi*wfs(wf).fc*wfs(wf).time_raw);
            accum(adc).data{accum_idx} = resample(double(accum(adc).data{accum_idx}), param.wfs(1).ft_dec(1), param.wfs(1).ft_dec(2));
          else
            accum(adc).data{accum_idx} = ifft(fft(accum(adc).data{accum_idx}) .* ifftshift(param.proc.ft_wind(length(accum(adc).data{accum_idx}))));
          end
          if ~param.load.wf_adc_comb.en
            % Regular wf-adc pair loading: no combining wf-adc pairs in fast-time
            if param.proc.combine_rx
              g_data{img_idx}(:,out_idx) = g_data{img_idx}(:,out_idx) + accum(adc).data{accum_idx} / param.proc.presums / size(param.load.imgs{img_idx},1);
            else
              g_data{img_idx}(:,out_idx,wf_adc_idx) = accum(adc).data{accum_idx} / param.proc.presums;
            end

          else
            % Combine wf-adc pairs in fast-time
            if accum(adc).img_comb_idx(accum_idx) == 1
              tmp2{wf_adc_idx} = zeros(param.load.wf_adc_comb.Nt_orig,1);
              tmp2{wf_adc_idx}(1:param.load.wf_adc_comb.rbins(1,out_idx)) ...
                = accum(adc).data{accum_idx}(1:param.load.wf_adc_comb.rbins(1,out_idx)) / param.proc.presums;
              %               g_data{img_idx}(1:param.load.wf_adc_comb.rbins(1,out_idx),out_idx,wf_adc_idx) ...
              %                 = accum(board+1).data{accum_idx}(1:param.load.wf_adc_comb.rbins(1,out_idx)) / param.proc.presums;
            elseif accum(adc).img_comb_idx(accum_idx) == 2
              tmp2{wf_adc_idx}(param.load.wf_adc_comb.rbins(1,out_idx)+1:end) ...
                = accum(adc).data{accum_idx}(param.load.wf_adc_comb.rbins(2,out_idx):end) / param.proc.presums;
              g_data{img_idx}(:,out_idx,wf_adc_idx) = tmp2{wf_adc_idx}(param.load.wf_adc_comb.keep_bins);
              %               g_data{img_idx}(param.load.wf_adc_comb.rbins(1,out_idx)+1:end,out_idx,wf_adc_idx) ...
              %                 = accum(board+1).data{accum_idx}(param.load.wf_adc_comb.rbins(2,out_idx):end) / param.proc.presums;
            end
          end
        end
      end
      if boundary_crossing
        rec=rec+1;
        boundary_crossing=false;
      end

    end
else %read incoherent data
    while rec < total_rec;
      % Get the filename
      fn_idx = param.load.file_idx{adc_idx}(rec);
      fn = param.load.filenames{adc_idx}{fn_idx};

      % Get number of records to load
      num_rec = find(param.load.file_idx{adc_idx}(rec:total_rec ) > fn_idx, 1) - 1;

      % Check to see if we crossed a file boundary
      if isempty(num_rec)
        % Not crossing a file boundary
        num_rec = length(rec:total_rec);
        boundary_crossing = false;
      else
        % Crossing a file boundary
        boundary_crossing = true;
      end

      % ===============================================================
      % Load the records
      % ===============================================================

      fprintf('Loading data from %s\n', fn);
      [fid,msg] = fopen(fn, 'r');
      if fid<0
        error('File open failed (%s)\n%s',fn, msg);
      end
      file_record_length=size(icards_get_data(fn,1),2);

      for rec = rec:rec+num_rec-1

        % Load data record from file
        fseek(fid, param.load.offset{1}(rec), 'bof');
        rec_data = [fread(fid, param.load.rec_data_size/sample_size, 'uint16')];%read coherent data      
        rec_data = double(rec_data);

        % ===============================================================
        % Process record
        % ===============================================================
        accum_idx = 1; %icards just has one waveform
        if num_accum == 0
          accum.data{accum_idx} = rec_data;
        else
          accum.data{accum_idx} = accum.data{accum_idx} + rec_data;
        end

        num_accum = num_accum + 1;
        if num_accum >= param.proc.presums
          num_accum = 0;
          out_idx = out_idx + 1;

          wf = accum(adc).wf(accum_idx);
          img_idx = accum(adc).img_idx(accum_idx);
          wf_adc_idx = accum(adc).wf_adc_idx(accum_idx);

          % Apply channel compensation
          if ~param.proc.raw_data
            chan_equal = 10.^(param.radar.wfs(wf).chan_equal_dB(param.radar.wfs(wf).rx_paths(adc))/20) ...
              .* exp(1i*param.radar.wfs(wf).chan_equal_deg(param.radar.wfs(wf).rx_paths(adc))/180*pi);
            accum(adc).data{accum_idx} = accum(adc).data{accum_idx}/chan_equal;
            accum(adc).data{accum_idx} = accum(adc).data{accum_idx}/wfs(wf).adc_gains(adc);
          end
          
          if param.proc.pulse_comp
            % ===========================================================
            % Do pulse compression
            % Apply blank (only should enable if sidelobe problems present)
            if ~isempty(wfs(wf).blank)
              % accum(adc).data{accum_idx}(wfs(wf).time_raw>wfs(wf).blank(1) & wfs(wf).time_raw<wfs(wf).blank(2)) = 0;
              accum(adc).data{accum_idx}(wfs(wf).time_raw-param.radar.wfs(wf).Tsys(adc) <= param.surface(rec) + wfs(wf).blank) = 0;
            end
            % Apply matched filter
            % Zero pad front: (the standard)
            accum(adc).data{accum_idx} = fft([zeros(wfs(wf).pad_length,1); accum(adc).data{accum_idx}], wfs(wf).Nt_pc);
            % Zero pad end: (debug only)
            %accum(adc).data{accum_idx} = fft(accum(adc).data{accum_idx}, wfs(wf).Nt_pc);
            accum(adc).data{accum_idx} = ifft(accum(adc).data{accum_idx}(wfs(wf).freq_inds) ...
              .* wfs(wf).ref{adc}(wfs(wf).freq_inds));
            if wfs(wf).dc_shift ~= 0
              % Correct for small frequency offset caused by selecting bins from
              % frequency domain as the method for down conversion
              accum(adc).data{accum_idx} = accum(adc).data{accum_idx}.*exp(-1i*2*pi*wfs(wf).dc_shift*wfs(wf).time);
            end
          elseif param.proc.ft_dec
            accum(adc).data{accum_idx} = fft(accum(adc).data{accum_idx},wfs(wf).Nt_raw);
            accum(adc).data{accum_idx} = ifft(accum(adc).data{accum_idx}(wfs(wf).freq_inds));
            if wfs(wf).dc_shift ~= 0
              % Correct for small frequency offset caused by selecting bins from
              % frequency domain as the method for down conversion
              accum(adc).data{accum_idx} = accum(adc).data{accum_idx}.*exp(-1i*2*pi*wfs(wf).dc_shift*wfs(wf).time);
            end
          end
          if ~param.load.wf_adc_comb.en
            % Regular wf-adc pair loading: no combining wf-adc pairs in fast-time
            if param.proc.combine_rx
              g_data{img_idx}(:,out_idx) = g_data{img_idx}(:,out_idx) + accum(adc).data{accum_idx} / param.proc.presums / size(param.load.imgs{img_idx},1);
            else
              g_data{img_idx}(:,out_idx,wf_adc_idx) = accum(adc).data{accum_idx} / param.proc.presums;
            end

          else
            % Combine wf-adc pairs in fast-time
            if accum(adc).img_comb_idx(accum_idx) == 1
              tmp2{wf_adc_idx} = zeros(param.load.wf_adc_comb.Nt_orig,1);
              tmp2{wf_adc_idx}(1:param.load.wf_adc_comb.rbins(1,out_idx)) ...
                = accum(adc).data{accum_idx}(1:param.load.wf_adc_comb.rbins(1,out_idx)) / param.proc.presums;
              %               g_data{img_idx}(1:param.load.wf_adc_comb.rbins(1,out_idx),out_idx,wf_adc_idx) ...
              %                 = accum(board+1).data{accum_idx}(1:param.load.wf_adc_comb.rbins(1,out_idx)) / param.proc.presums;
            elseif accum(adc).img_comb_idx(accum_idx) == 2
              tmp2{wf_adc_idx}(param.load.wf_adc_comb.rbins(1,out_idx)+1:end) ...
                = accum(adc).data{accum_idx}(param.load.wf_adc_comb.rbins(2,out_idx):end) / param.proc.presums;
              g_data{img_idx}(:,out_idx,wf_adc_idx) = tmp2{wf_adc_idx}(param.load.wf_adc_comb.keep_bins);
              %               g_data{img_idx}(param.load.wf_adc_comb.rbins(1,out_idx)+1:end,out_idx,wf_adc_idx) ...
              %                 = accum(board+1).data{accum_idx}(param.load.wf_adc_comb.rbins(2,out_idx):end) / param.proc.presums;
            end
          end
        end
      end
      if boundary_crossing
        rec=rec+1;
        boundary_crossing=false;
      end

    end
    
end
fclose(fid);
% one-sided tukey window
[ window_matrix ] = onesided_window( whole_param.radar.icards.td_window,whole_param.radar.icards.td_window_side,size(g_data{1}) );%using one-sided window to avoid too big value at the start or end of a record
g_data{img_idx}=g_data{img_idx}.*window_matrix;
% if strcmpi(whole_param.radar.icards.data_type,'coherent') %incoherent data does not need burst noise detection or slow time noise removal
    % burst noise detection
%     [g_data{img_idx}]=icards_burst_noise_detection(g_data{img_idx});
    % DC offset removals
    g_data{img_idx} = g_data{img_idx} - repmat(fir_dec(mean(g_data{img_idx}(round(size(g_data{1},1)*whole_param.radar.icards.noise_cal_rng):end,:),1), whole_param.radar.icards.DC_filter), [size(g_data{img_idx},1) 1]);
% end
return;

% ==============================================================
% accum = build_img_load_struct(imgs, adcs)
% ==============================================================
function accum = build_img_load_struct(imgs, adcs)

for img_idx = 1:length(imgs)
  for wf_adc_idx = 1:size(imgs{img_idx},1)
    if ~ismember(imgs{img_idx}(wf_adc_idx,2), adcs)
      error('ADC %d is invalid', imgs{img_idx}(wf_adc_idx,2));
    end
  end
end

% Build accum structure
for adc = adcs
  accum(adc).wf = [];
  accum(adc).wf_adc_idx = [];
  accum(adc).img_idx = [];
  accum(adc).img_comb_idx = [];
  for img_idx = 1:length(imgs)
    for wf_adc_idx = 1:size(imgs{img_idx},1)
      for adc_column = 2:2:size(imgs{img_idx},2)
        if imgs{img_idx}(wf_adc_idx,2) == adc
          accum(adc).wf(end+1) = imgs{img_idx}(wf_adc_idx,adc_column-1);
          accum(adc).wf_adc_idx(end+1) = wf_adc_idx;
          accum(adc).img_idx(end+1) = img_idx;
          accum(adc).img_comb_idx(end+1) = adc_column/2;
        end
      end
    end
  end
end

return;


