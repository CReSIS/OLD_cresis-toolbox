function [hdr,data] = data_load(param,records,states)
% [hdr,data] = data_load(param,records,states)
%
% https://ops.cresis.ku.edu/wiki/index.php/Data_load#data_load.m
%
% Author: John Paden

wfs = param.radar.wfs;

%% Preallocate data
% ===================================================================
total_rec = param.load.recs(end)-param.load.recs(1)+1;
Nx = floor(total_rec/param.load.presums);
data = cell(size(param.load.imgs));
hdr = [];
hdr.bad_rec = cell(size(param.load.imgs));
for img = 1:length(param.load.imgs)
  wf = abs(param.load.imgs{img}(1,1));
  Nt = wfs(wf).Nt_raw;
  Nc = size(param.load.imgs{img},1);
  % data{img}: may eventually be complex(zeros(Nt,Nx,Nc,'single')) after
  % the first record/range-line is written
  data{img} = zeros(0,Nx,Nc,'single');
  hdr.bad_rec{img} = zeros(1,Nx,Nc,'uint8');
  hdr.nyquist_zone_hw{img} = zeros(1,Nx,'uint8');
  hdr.nyquist_zone_signal{img} = nan(1,Nx);
  hdr.DDC_dec{img} = ones(1,Nx,'double');
  hdr.DDC_freq{img} = zeros(1,Nx,'double');
  hdr.Nt{img} = zeros(1,Nx,'double');
  hdr.t0_raw{img} = zeros(1,Nx,'double');
  hdr.t_ref{img} = zeros(1,Nx,'double');
end
nyquist_zone_hw = zeros(1,param.load.presums);
nyquist_zone_signal = NaN;
DDC_dec = ones(1,param.load.presums);
DDC_freq = zeros(1,param.load.presums);
Nt = zeros(1,param.load.presums);
t0 = zeros(1,param.load.presums);
t_ref = zeros(1,param.load.presums);

%% Endian mode
% ===================================================================
[~,~,endian] = computer;
if any(param.records.file.version==[9 10 103 412])
  if endian == 'B'
    % IEEE little endian files with IEEE big endian computer
    swap_bytes_en = true;
  else
    swap_bytes_en = false;
  end
else
  if endian == 'L'
    % IEEE big endian files with IEEE little endian computer
    swap_bytes_en = true;
  else
    swap_bytes_en = false;
  end
end

%% Waveform Map (t_ref calculation)
% ===================================================================
if any(param.records.file.version==[8])
  % Handles bit shift error in waveform_ID field
  waveform_ID_map = [typecast(uint8('     000'),'uint64')
    typecast(uint8('     090'),'uint64')
    typecast(uint8('     120'),'uint64')
    typecast(uint8('     330'),'uint64')
    typecast(uint8('     000')*2,'uint64')
    typecast(uint8('     090')*2,'uint64')
    typecast(uint8('     120')*2,'uint64')
    typecast(uint8('     330')*2,'uint64')
    typecast(uint8('OIB_FMCW'),'uint64')
    typecast(uint8('OIB_FMCW')*2,'uint64')
    ];
  waveform_ID_t_ref = [0e-6 9e-6 12e-6 33e-6 0e-6 9e-6 12e-6 33e-6 0e-6 0e-6];
  
end

%% Load data
% ===================================================================
for state_idx = 1:length(states)
  state = states(state_idx);
  file_data_last_file = [];
  board = state.board;
  board_idx = state.board_idx;
  fid = 0;
  out_rec = 0;
  num_accum = 0;
  num_presum_records = 0;

  % Find the file index for each record. If records.offset has extra
  % entries, get file_idxs for these too.
  file_idxs = relative_rec_num_to_file_idx_vector( ...
    param.load.recs+[0 size(records.offset,2)-(1+diff(param.load.recs))],records.relative_rec_num{board_idx});
  
  rec = 1;
  while rec <= total_rec
    
    %% Load in a file
    if records.offset(board_idx,rec) ~= -2^31
      % Determine which file has the current record
      file_idx = file_idxs(rec);
      if ~isempty(file_data_last_file)
        % Part of current record has been loaded from previous file
      elseif records.offset(board_idx,rec) < 0 && records.offset(board_idx,rec) ~= -2^31
        % Record offset is negative, but not -2^31: this means the record
        % started in the previous file.
        file_idx = file_idx - 1;
      end
      
      % Get the file's name
      adc = state.adc(1);
      fn_name = records.relative_filename{board_idx}{file_idx};
      [fn_dir] = get_segment_file_list(param,board_idx);
      fn = fullfile(fn_dir,fn_name);

      % Open the file
      fprintf('  Open %s (%s)\n', fn, datestr(now));
      [fid,msg] = fopen(fn, 'rb');
      if fid <= 0
        error('File open failed (%s)\n%s',fn, msg);
      end

      % Seek to the current record's position in the file
      if ~isempty(file_data_last_file)
        file_data_offset = records.offset(board_idx,rec);
      elseif records.offset(board_idx,rec) < 0
        finfo = dir(fn);
        file_data_offset = finfo.bytes + records.offset(board_idx,rec);
        fseek(fid,file_data_offset,-1);
      else
        file_data_offset = records.offset(board_idx,rec);
        fseek(fid,file_data_offset,-1);
      end
      
      % Load the rest of the file into memory
      file_data = [file_data_last_file(:); fread(fid,inf,'uint8=>uint8')];
      fclose(fid);
    else
      file_idx = -1;
    end

    %% Pull out records from this file
    while rec <= total_rec
      if records.offset(board_idx,rec) ~= -2^31
        if records.offset(board_idx,rec) < 0
          if isempty(file_data_last_file)
            % Record offset is negative and so represents a record that
            % bridges into the next file
            file_data_last_file = file_data(end+records.offset(board_idx,rec)+1:end);
            break
          else
            file_data_last_file = [];
          end
        end
        if file_idxs(rec) > file_idx
          break;
        end
        
        % Extract next record (determine its relative position in the
        % file_data memory block
        rec_offset = records.offset(board_idx,rec) - file_data_offset;

        % record_mode==1: Find the size of the current record
        if rec < length(file_idxs) && file_idxs(rec+1) == file_idxs(rec)
          % Next record is in this file, rec size is set to the offset to
          % the next record
          rec_size = records.offset(board_idx,rec+1) - records.offset(board_idx,rec);
        else
          % Next record is in the next file, rec_size is set to the rest of
          % the data in this file
          rec_size = (length(file_data)-file_data_offset) - records.offset(board_idx,rec);
        end
        
        % Process all adc-wf pairs in this record
        missed_wf_adc = false;
        for accum_idx = 1:length(state.wf)
          adc = state.adc(accum_idx);
          wf = state.wf(accum_idx);
          mode_latch = state.mode(accum_idx);
          subchannel = state.subchannel(accum_idx);
          
          % Read in headers for this waveform
          % ---------------------------------------------------------------
          
          if wfs(wf).quantization_to_V_dynamic
            if param.records.file.version == 407
              bit_shifts = -typecast(file_data(rec_offset + wfs(wf).offset - 4),'int8');
              % Apply dynamic bit shifts
              quantization_to_V_adjustment = 2^(bit_shifts - wfs(wf).bit_shifts(adc));
            end
          else
            quantization_to_V_adjustment = 1;
          end

          % Read in headers for this record
          if any(param.records.file.version == [3 5 7 8])
            % Number of fast-time samples Nt, and start time t0
            if swap_bytes_en
              start_idx = 2*double(swapbytes(typecast(file_data(rec_offset+37:rec_offset+38), 'uint16'))) + wfs(wf).time_raw_trim(1);
              stop_idx = 2*double(swapbytes(typecast(file_data(rec_offset+39:rec_offset+40), 'uint16'))) - wfs(wf).time_raw_trim(2);
            else
              start_idx = 2*double(typecast(file_data(rec_offset+37:rec_offset+38), 'uint16')) + wfs(wf).time_raw_trim(1);
              stop_idx = 2*double(typecast(file_data(rec_offset+39:rec_offset+40), 'uint16')) - wfs(wf).time_raw_trim(2);
            end
            if param.records.file.version == 8
              Nt(num_accum+1) = stop_idx - start_idx;
              wfs(wf).Nt_raw = Nt(num_accum+1);
            else
              % NCO frequency
              if swap_bytes_en
                DDC_freq(num_accum+1) = double(-swapbytes(typecast(file_data(rec_offset+43:rec_offset+44),'uint16')));
              else
                DDC_freq(num_accum+1) = double(-typecast(file_data(rec_offset+43:rec_offset+44),'uint16'));
              end
              if param.records.file.version == 3
                DDC_freq(num_accum+1) = DDC_freq(num_accum+1) / 2^15 * wfs(wf).fs_raw * 2 - 62.5e6;
              elseif any(param.records.file.version == [5 7])
                DDC_freq(num_accum+1) = DDC_freq(num_accum+1) / 2^15 * wfs(wf).fs_raw * 2;
              end
              
              DDC_dec(num_accum+1) = 2^(double(file_data(rec_offset+46))+1);
              raw_or_DDC = file_data(rec_offset + wfs(wf).offset +48);
              if raw_or_DDC
                Nt(num_accum+1) = (stop_idx - start_idx);
              else
                Nt(num_accum+1) = floor((stop_idx - start_idx) / DDC_dec(num_accum+1));
              end
              wfs(wf).Nt_raw = Nt(num_accum+1);
            end
            t0(num_accum+1) = start_idx/wfs(wf).fs_raw;
            
            if param.records.file.version == 8
              % Debug: char(file_data(rec_offset+41:rec_offset+48).')
              % No swapbytes should be necessary for this typecast because
              % it is actually a string of 8 characters.
              waveform_ID = typecast(file_data(rec_offset+41:rec_offset+48), 'uint64');
              waveform_ID_map_idx = find(waveform_ID_map == waveform_ID,1);
              if isempty(waveform_ID_map_idx)
                error('%ld waveform_ID not found in waveform_ID_map.',waveform_ID);
              end
              t_ref(num_accum+1) = wfs(wf).t_ref + waveform_ID_t_ref(waveform_ID_map_idx);
            else
              % Reference deramp time delay, t_ref
              t_ref(num_accum+1) = wfs(wf).t_ref;
            end
            
            % Bit shifts
            if wfs(wf).quantization_to_V_dynamic
              bit_shifts = double(-typecast(file_data(rec_offset+36),'int8'));
              quantization_to_V_adjustment = 2^(bit_shifts - wfs(wf).bit_shifts(adc));
            end
            
            % Nyquist zone
            if param.records.file.version == 8
              nyquist_zone_hw(num_accum+1) = file_data(rec_offset+34);
            elseif any(param.records.file.version == [3 5 7])
              nyquist_zone_hw(num_accum+1) = file_data(rec_offset+45);
            end
          end
          nyquist_zone_signal = nyquist_zone_hw(1);
          if isfield(records.settings,'nyquist_zone') && ~isnan(records.settings.nyquist_zone(rec))
            nyquist_zone_signal = records.settings.nyquist_zone(rec);
          end
          if size(data{state.img(accum_idx)},1) < wfs(wf).Nt_raw
            % Force data output to grow to the current record size
            data{state.img(accum_idx)}(wfs(wf).Nt_raw,1,1) = 0;
          end
          

          % Extract waveform for this wf-adc pair
          switch wfs(wf).record_mode
            case 0
              % Read in standard fixed record
              %  - Supports interleaved IQ samples
              %  - Supports arbitrary sample types
              %  - Supports interleaved data channels ("adcs")
              start_bin = 1+rec_offset + wfs(wf).offset + wfs(wf).time_raw_trim(1)*wfs(wf).adc_per_board*wfs(wf).sample_size;
              stop_bin = start_bin + wfs(wf).Nt_raw*wfs(wf).adc_per_board*wfs(wf).sample_size-1;
              if swap_bytes_en
                tmp = single(swapbytes(typecast(file_data(start_bin : stop_bin), wfs(wf).sample_type)));
              else
                tmp = single(typecast(file_data(start_bin : stop_bin), wfs(wf).sample_type));
              end
              if wfs(wf).complex
                if wfs(wf).conjugate
                  tmp = tmp(1:2:end) - 1i*tmp(2:2:end);
                else
                  tmp = tmp(1:2:end) + 1i*tmp(2:2:end);
                end
              end
              adc_offset = mod(adc-1,wfs(wf).adc_per_board);
              tmp = tmp(1+adc_offset : wfs(wf).adc_per_board : end);
              if param.records.file.version ~= 408
                tmp_data{adc,wf} = tmp;
              else
                % 8 sample interleave for file_version 408
                tmp_data{adc,wf}(1:8:length(tmp)) = tmp(1:8:end);
                tmp_data{adc,wf}(5:8:length(tmp)) = tmp(5:8:end);
                tmp_data{adc,wf}(2:8:length(tmp)) = tmp(2:8:end);
                tmp_data{adc,wf}(6:8:length(tmp)) = tmp(6:8:end);
                tmp_data{adc,wf}(3:8:length(tmp)) = tmp(3:8:end);
                tmp_data{adc,wf}(7:8:length(tmp)) = tmp(7:8:end);
                tmp_data{adc,wf}(4:8:length(tmp)) = tmp(4:8:end);
                tmp_data{adc,wf}(8:8:length(tmp)) = tmp(8:8:end);
              end
              
            case 1
              % Read in RSS dynamic record
              sub_rec_offset = 0;
              missed_wf_adc = true;
              while sub_rec_offset < rec_size
                total_offset = rec_offset + sub_rec_offset;
                if swap_bytes_en
                  radar_header_type = mod(swapbytes(typecast(file_data(total_offset+(9:12)),'uint32')),2^31); % Ignore MSB
                  radar_header_len = double(swapbytes(typecast(file_data(total_offset+(13:16)),'uint32')));
                  radar_profile_length = double(swapbytes(typecast(file_data(total_offset+radar_header_len+(21:24)),'uint32')));
                else
                  radar_header_type = mod(typecast(file_data(total_offset+(9:12)),'uint32'),2^31); % Ignore MSB
                  radar_header_len = double(typecast(file_data(total_offset+(13:16)),'uint32'));
                  radar_profile_length = double(typecast(file_data(total_offset+radar_header_len+(21:24)),'uint32'));
                end
                if any(radar_header_type == [5 16 23])
                  if mode_latch == typecast(file_data(total_offset+17),'uint8') ...
                      && subchannel == typecast(file_data(total_offset+18),'uint8')
                    % This matches the mode and subchannel that we need
                    is_IQ = 0;
                    if length(file_data) < total_offset+24+radar_header_len+radar_profile_length
                      % Unexpected end of file, so we missed the record
                      missed_wf_adc = true;
                      break;
                    end
                    if swap_bytes_en
                      radar_profile_format = swapbytes(typecast(file_data(total_offset+radar_header_len+(17:20)),'uint32'));
                      switch radar_profile_format
                        case 0 % 0x00000
                          tmp_data{adc,wf} = single(swapbytes(typecast(file_data(total_offset+24+radar_header_len+(1:radar_profile_length)),'int16')));
                        case 65536 % 0x10000
                          tmp_data{adc,wf} = single(swapbytes(typecast(file_data(total_offset+24+radar_header_len+(1:radar_profile_length)),'int16')));
                          is_IQ = 1;
                        case 131072 % 0x20000
                          tmp_data{adc,wf} = single(swapbytes(typecast(file_data(total_offset+24+radar_header_len+(1:radar_profile_length)),'int32')));
                          is_IQ = 1;
                        case 196608 % 0x30000
                          tmp_data{adc,wf} = single(swapbytes(typecast(file_data(total_offset+24+radar_header_len+(1:radar_profile_length)),'single')));
                          is_IQ = 1;
                      end
                    else
                      radar_profile_format = typecast(file_data(total_offset+radar_header_len+(17:20)),'uint32');
                      switch radar_profile_format
                        case 0 % 0x00000
                          tmp_data{adc,wf} = single(typecast(file_data(total_offset+24+radar_header_len+(1:radar_profile_length)),'int16'));
                        case 65536 % 0x10000
                          tmp_data{adc,wf} = single(typecast(file_data(total_offset+24+radar_header_len+(1:radar_profile_length)),'int16'));
                          is_IQ = 1;
                        case 131072 % 0x20000
                          tmp_data{adc,wf} = single(typecast(file_data(total_offset+24+radar_header_len+(1:radar_profile_length)),'int32'));
                          is_IQ = 1;
                        case 196608 % 0x30000
                          tmp_data{adc,wf} = single(typecast(file_data(total_offset+24+radar_header_len+(1:radar_profile_length)),'single'));
                          is_IQ = 1;
                      end
                    end
                    if is_IQ
                      if wfs(wf).conjugate
                        tmp_data{adc,wf} = tmp_data{adc,wf}(1:2:end) - 1i*tmp_data{adc,wf}(2:2:end);
                      else
                        tmp_data{adc,wf} = tmp_data{adc,wf}(1:2:end) + 1i*tmp_data{adc,wf}(2:2:end);
                      end
                    end
                    % Found the record so break
                    missed_wf_adc = false;
                    break;
                  end
                else
                  % Unsupported radar header type, skip this record
                  break;
                end
                sub_rec_offset = sub_rec_offset + 24 + radar_header_len + radar_profile_length;
              end
          end
          
          if missed_wf_adc
            break;
          end
          
          if quantization_to_V_adjustment ~= 1 && ~param.load.raw_data
            % Convert from quantization to voltage at the receiver input for the
            % maximum gain case:
            %  1. fast time gains less than the maximum for this record will be
            %     compensated for in the next step
            %  2. antenna effects not considered at this step
            tmp_data{adc,wf} = tmp_data{adc,wf} * quantization_to_V_adjustment;
            
          end
          
          % Accumulate (presum)
          if num_accum == 0
            state.data{accum_idx} = tmp_data{adc,wf};
          else
            state.data{accum_idx} = state.data{accum_idx} + tmp_data{adc,wf};
          end
        end
        if ~missed_wf_adc
          num_accum = num_accum + 1;
        end
      end
      
      % Store to output if number of presums is met
      num_presum_records = num_presum_records + 1;
      if num_presum_records >= param.load.presums
        out_rec = out_rec + 1;
        for accum_idx = 1:length(state.wf)
          % Sum up wf-adc sum pairs until done
          if num_accum >= 1
            switch state.wf_adc_sum_cmd(accum_idx)
              case 0
                state.data{accum_idx} = state.wf_adc_sum(accum_idx)*state.data{accum_idx};
                continue;
              case 1
                state.data{accum_idx} = state.data{accum_idx} ...
                  + state.wf_adc_sum(accum_idx)*state.data{accum_idx};
                continue;
              case 2
                state.data{accum_idx} = state.data{accum_idx} ...
                  + state.wf_adc_sum(accum_idx)*state.data{accum_idx};
              case 3
                state.data{accum_idx} = state.wf_adc_sum(accum_idx)*state.data{accum_idx};
            end
          end

          % Store to output
          if num_accum < num_presum_records*wfs(wf).presum_threshold ...
              || any(nyquist_zone_hw ~= nyquist_zone_hw(1)) ...
              || any(DDC_dec ~= DDC_dec(1)) ...
              || any(DDC_freq ~= DDC_freq(1)) ...
              || any(Nt ~= Nt(1)) ...
              || any(t0 ~= t0(1)) ...
              || any(t_ref ~= t_ref(1))
            % Too few presums, mark as bad record
            % Or a parameter changed within the presum block
            data{state.img(accum_idx)}(:,out_rec,state.wf_adc(accum_idx)) = 0;
            hdr.bad_rec{state.img(accum_idx)}(1,out_rec,state.wf_adc(accum_idx)) = 1;
          else
            data{state.img(accum_idx)}(:,out_rec,state.wf_adc(accum_idx)) = state.data{accum_idx};
          
            hdr.nyquist_zone_hw{img}(out_rec) = nyquist_zone_hw(1);
            hdr.nyquist_zone_signal{img}(out_rec) = nyquist_zone_signal;
            hdr.DDC_dec{img}(out_rec) = DDC_dec(1);
            hdr.DDC_freq{img}(out_rec) = DDC_freq(1);
            hdr.Nt{img}(out_rec) = Nt(1);
            hdr.t0_raw{img}(out_rec) = t0(1);
            hdr.t_ref{img}(out_rec) = t_ref(1);
            
            if any(isfinite(wfs(wf).blank))
              % Blank data
              %  - Blank is larger of two numbers passed in through radar worksheet blank parameter:
              %   Number 1 is added to surface time delay and is usually equal to pulse duration
              %   Number 2 is unmodified and is usually equal to hardware blank setting
              %   Set either number to -inf to disable
              blank_time = max(records.surface(rec) + wfs(wf).blank(1),wfs(wf).blank(2));
              data{state.img(accum_idx)}(wfs(wf).time_raw-param.radar.wfs(wf).Tsys(adc) <= blank_time,out_rec,state.wf_adc(accum_idx)) = 0;
            end
          end
        end
        
        % Reset counters
        num_presum_records = 0;
        num_accum = 0;
      end
      % Increment record counter
      rec = rec + 1;
    end
  end
  
end

if ~param.load.raw_data
  for img = 1:length(param.load.imgs)
    for wf_adc = 1:size(data{img},3)
      wf = param.load.imgs{img}(wf_adc,1);
      adc = param.load.imgs{img}(wf_adc,2);
      
      % Apply channel compensation, presum normalization, and constant
      % receiver gain compensation
      chan_equal = 10.^(param.radar.wfs(wf).chan_equal_dB(param.radar.wfs(wf).rx_paths(adc))/20) ...
        .* exp(1i*param.radar.wfs(wf).chan_equal_deg(param.radar.wfs(wf).rx_paths(adc))/180*pi);
      mult_factor = single(wfs(wf).quantization_to_V(adc)/param.load.presums/10.^(wfs(wf).adc_gains_dB(adc)/20)/chan_equal);
      data{img}(:,:,wf_adc) = mult_factor * data{img}(:,:,wf_adc);
      
      % Compensate for receiver gain applied before ADC quantized the signal
      %  - For time varying receiver gain, the convention is to compensate
      %    to the maximum receiver gain and use the wfs(wf).gain parameter
      %    to vary the gain relative to that.
      % Apply fast-time varying gain if defined
      if ~isempty(wfs(wf).gain)
        data{img}(:,:,wf_adc) = bsxfun(@times,data{img}(:,:,wf_adc),interp1(wfs(wf).gain.Time, wfs(wf).gain.Gain, wfs(wf).time_raw(1:wfs(wf).Nt_raw)));
      end
      
      % Apply time varying channel compensation
      if ~isempty(wfs(wf).chan_equal)
        cdf_fn_dir = fileparts(ct_filename_out(param,wfs(wf).chan_equal, ''));
        cdf_fn = fullfile(cdf_fn_dir,sprintf('chan_equal_%s_wf_%d_adc_%d.nc', param.day_seg, wf, adc));
        
        finfo = ncinfo(cdf_fn);
        % Determine number of records and set recs(1) to this
        Nt = finfo.Variables(find(strcmp('chan_equal',{finfo.Variables.Name}))).Size(2);
        
        chan_equal = [];
        chan_equal.gps_time = ncread(cdf_fn,'gps_time');
        recs = find(chan_equal.gps_time > records.gps_time(1) - 100 & chan_equal.gps_time < records.gps_time(end) + 100);
        chan_equal.gps_time = chan_equal.gps_time(recs);
        
        chan_equal.chan_equal = ncread(cdf_fn,'chan_equalI',[recs(1) 1],[recs(end)-recs(1)+1 Nt]) ...
          + 1i*ncread(cdf_fn,'chan_equalQ',[recs(1) 1],[recs(end)-recs(1)+1 Nt]);
        
        data{img}(1:wfs(wf).Nt_raw,:,wf_adc) = data{img}(1:wfs(wf).Nt_raw,:,wf_adc) ...
          .*interp1(reshape(chan_equal.gps_time,[numel(chan_equal.gps_time) 1]),chan_equal.chan_equal,records.gps_time,'linear','extrap').';
      end
    end
  end
end

%% Add record information
% =========================================================================
hdr.gps_source = records.gps_source;
hdr.gps_time = fir_dec(records.gps_time,param.load.presums);
hdr.surface = fir_dec(records.surface,param.load.presums);
  
%% Create trajectories

% Create reference trajectory (rx_path == 0, tx_weights = [])
trajectory_param = struct('gps_source',records.gps_source, ...
  'season_name',param.season_name,'radar_name',param.radar_name,'rx_path', 0, ...
  'tx_weights', [], 'lever_arm_fh', param.radar.lever_arm_fh);
hdr.ref = trajectory_with_leverarm(records,trajectory_param);

for img = 1:length(param.load.imgs)
  
  for wf_adc = 1:size(param.load.imgs{img},1)
    wf = param.load.imgs{img}(wf_adc,1);
    adc = param.load.imgs{img}(wf_adc,2);
    
    if isempty(param.radar.lever_arm_fh)
      hdr.records{img,wf_adc} = records;
      hdr.records{img,wf_adc}.lat = fir_dec(records.lat,param.load.presums);
      hdr.records{img,wf_adc}.lon = fir_dec(records.lon,param.load.presums);
      hdr.records{img,wf_adc}.elev = fir_dec(records.elev,param.load.presums);
      hdr.records{img,wf_adc}.roll = fir_dec(records.roll,param.load.presums);
      hdr.records{img,wf_adc}.pitch = fir_dec(records.pitch,param.load.presums);
      hdr.records{img,wf_adc}.heading = fir_dec(records.heading,param.load.presums);
    else
      % Create actual trajectory
      trajectory_param = struct('gps_source',records.gps_source, ...
        'season_name',param.season_name,'radar_name',param.radar_name, ...
        'rx_path', wfs(wf).rx_paths(adc), ...
        'tx_weights', wfs(wf).tx_weights, 'lever_arm_fh', param.radar.lever_arm_fh);
      hdr.records{img,wf_adc} = trajectory_with_leverarm(records,trajectory_param);
      hdr.records{img,wf_adc}.lat = fir_dec(hdr.records{img,wf_adc}.lat,param.load.presums);
      hdr.records{img,wf_adc}.lon = fir_dec(hdr.records{img,wf_adc}.lon,param.load.presums);
      hdr.records{img,wf_adc}.elev = fir_dec(hdr.records{img,wf_adc}.elev,param.load.presums);
      hdr.records{img,wf_adc}.roll = fir_dec(hdr.records{img,wf_adc}.roll,param.load.presums);
      hdr.records{img,wf_adc}.pitch = fir_dec(hdr.records{img,wf_adc}.pitch,param.load.presums);
      hdr.records{img,wf_adc}.heading = fir_dec(hdr.records{img,wf_adc}.heading,param.load.presums);
    end
    hdr.records{img,wf_adc}.offset = hdr.records{img,wf_adc}.offset(:,1:param.load.presums:end-param.load.presums+1);
    hdr.records{img,wf_adc} = rmfield(hdr.records{img,wf_adc},{'gps_time','surface'});
  end
end

