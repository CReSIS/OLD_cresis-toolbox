function [hdr,data] = data_load(param,records,wfs,states)
% [hdr,data] = data_load(param,records,wfs,states)

%% Preallocate data
% ===================================================================
total_rec = param.load.recs(end)-param.load.recs(1)+1;
Nx = floor(total_rec/param.load.presums);
data = cell(size(param.load.imgs));
hdr = cell(size(param.load.imgs));
for img = 1:length(param.load.imgs)
  wf = abs(param.load.imgs{img}(1,1));
  Nt = wfs(wf).Nt_raw;
  Nc = size(param.load.imgs{img},1);
  data{img} = zeros(Nt,Nx,Nc,'single');
  hdr{img}.bad_rec = zeros(Nx,Nc,'uint8');
end

%% Endian mode
% ===================================================================
if any(param.records.file_version==[9 411 412])
  file_mode = 'ieee-le';
else
  file_mode = 'ieee-be';
end


%% Load data
% ===================================================================
for state_idx = 1:length(states)
  state = states(state_idx);
  file_data_state = 0;
  file_data_last_file = [];
  board = state.board(state_idx);
  board_idx = state.board_idx(state_idx);
  fid = 0;
  out_rec = 0;
  num_accum = 0;
  num_presum_records = 0;
  
  file_idxs = relative_rec_num_to_file_idx_vector( ...
    param.load.recs,records.relative_rec_num{board_idx});
  
  rec = 1;
  while rec <= total_rec
    
    %% Load in a file
    if records.offset(board_idx,rec) ~= -2^31
      % Determine which file has the current record
      file_idx = file_idxs(rec);
      if records.offset(board_idx,rec) < 0 && records.offset(board_idx,rec) ~= -2^31
        % Record offset is negative, but not -2^31: this means the record
        % started in the previous file.
        file_idx = file_idx - 1;
      end
      
      % Get the file's name
      adc = state.adc(1);
      fn_name = records.relative_filename{board_idx}{file_idx};
      [fn_dir] = get_segment_file_list(param,adc);
      fn = fullfile(fn_dir,fn_name);

      % Open the file
      [fid,msg] = fopen(fn, 'r',file_mode);
      if fid <= 0
        error('File open failed (%s)\n%s',fn, msg);
      end

      % Seek to the current record's position in the file
      if records.offset(board_idx,rec) < 0
        if num_bytes==inf
          finfo = dir(fn)
        end
        file_data_offset = finfo.bytes + records.offset(board_idx,rec);
      else
        file_data_offset = records.offset(board_idx,rec);
      end
      fseek(fid,file_data_offset,-1);
      
      % Load the rest of the file into memory
      file_data = [file_data_last_file fread(fid,inf,'uint8=>uint8')];
    end
    file_data_state = 1;

    %% Pull out records from this file
    while file_data_state == 1 && rec <= total_rec
      
      if records.offset(board_idx,rec) ~= -2^31
        % Read in headers for this record
        % TO DO
        
        % Extract next record (determine its relative position in the
        % file_data memory block
        rec_offset = records.offset(board_idx,rec) - file_data_offset;
        
        % Check to see if enough data remains in the buffer for a new
        % record
        if size(records.offset,2) > rec
          rec_size = records.offset(board_idx,rec+1) - records.offset(board_idx,rec);
        elseif wfs(wf).record_mode ~= 1
          rec_size = wfs(wf).rec_size;
        end
        if numel(file_data) < rec_offset + rec_size
          % Not enough data left in this file to load the whole record,
          % go to read in next file state and store off the remainder of
          % this file's data
          file_data_last_file = file_data(rec_offset:end);
          file_data_state = 0;
          continue;
        end
        
        % Process all adc-wf pairs in this record
        for accum_idx = 1:length(state.wf)
          adc = state.adc(accum_idx);
          wf = state.wf(accum_idx);
          
          % Read in headers for this waveform
          % TO DO
          quantization_to_V_adjustment = 1;

          % Extract waveform for this wf-adc pair
          switch wfs(wf).record_mode
            case 0
              % Read in standard fixed record
              %  - Supports interleaved IQ samples
              %  - Supports arbitrary sample types
              %  - Supports interleaved data channels ("adcs")
              tmp = single(typecast(file_data(rec_offset + wfs(wf).offset + (0:wfs(wf).Nt_raw*wfs(wf).adc_per_board*wfs(wf).sample_size-1)), wfs(wf).sample_type));
              if wfs(wf).complex
                if wfs(wf).conjugate
                  tmp = tmp(1:2:end) + 1i*tmp(2:2:end);
                else
                  tmp = tmp(1:2:end) - 1i*tmp(2:2:end);
                end
              end
              adc_offset = mod(adc-1,wfs(wf).adc_per_board);
              tmp = tmp(1+adc_offset : wfs(wf).adc_per_board : end);
              tmp_data{adc,wf} = tmp;
              
            case 1
              % Read in RSS dynamic record
              
            case 2
              % Fixed records, 8 sample interleave for file_version 408
              tmp = single(typecast(file_data(rec_offset + wfs(wf).offset + (0:wfs(wf).Nt_raw*wfs(wf).sample_size-1)), wfs(wf).sample_type));
              tmp(1:8:end) = tmp(1:8:end);
              tmp(2:8:end) = tmp(5:8:end);
              tmp(3:8:end) = tmp(2:8:end);
              tmp(4:8:end) = tmp(6:8:end);
              tmp(5:8:end) = tmp(3:8:end);
              tmp(6:8:end) = tmp(7:8:end);
              tmp(7:8:end) = tmp(4:8:end);
              tmp(8:8:end) = tmp(8:8:end);
              tmp_data{adc,wf} = tmp;
              
          end
          
          if ~param.load.raw_data
            % Convert from quantization to voltage at the receiver input for the
            % maximum gain case:
            %  1. fast time gains less than the maximum for this record will be
            %     compensated for in the next step
            %  2. antenna effects not considered at this step
            tmp_data{adc,wf} = tmp_data{adc,wf} * wfs(wf).quantization_to_V * quantization_to_V_adjustment;
            
            % Apply fast-time varying gain
            if ~isempty(wfs(wf).gain)
              tmp_data{adc,wf} = tmp_data{adc,wf} .* interp1(wfs(wf).gain.Time, wfs(wf).gain.Gain, wfs(wf).time_raw(1:wfs(wf).Nt_raw));
            end
          end
          
          % Accumulate (presum)
          if num_accum == 0
            state.data{accum_idx} = tmp_data{adc,wf};
          else
            state.data{accum_idx} = state.data{accum_idx} + tmp_data{adc,wf};
          end
        end
        num_accum = num_accum + 1;
      end
      
      % Store to output if number of presums is met
      num_presum_records = num_presum_records + 1;
      if num_presum_records >= param.load.presums
        out_rec = out_rec + 1;
        for accum_idx = 1:length(state.wf)
          % Sum up wf-adc sum pairs until done
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
          
          % Store to output
          if num_accum < num_presum_records*wfs(wf).presum_threshold
            % Too few presums, mark as bad record
            data{state.img}(:,out_rec,state.wf_adc_idx(accum_idx)) = 0;
            hdr{state.img}.bad_rec(out_rec,state.wf_adc_idx(accum_idx)) = 1;
          else
            data{state.img(accum_idx)}(:,out_rec,state.wf_adc_idx(accum_idx)) = state.data{accum_idx} / num_accum;
          end
        end
        
        % Reset counters and increment record counter
        num_presum_records = 0;
        num_accum = 0;
        rec = rec + 1;
      end
    end
  end
  
end

return
