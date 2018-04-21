function [hdr,data] = data_load(param,wfs,records)
% [hdr,data] = data_load(param,records)

%% Preallocate data
% ===================================================================
total_rec = param.load.recs(end)-param.load.recs(1)+1;
Nx = floor(total_rec/param.proc.presums);
data = cell(size(param.load.imgs));
for img = 1:length(param.load.imgs)
  wf = abs(param.load.imgs{img}(1,1));
  if param.load.pulse_comp
    Nt = param.wfs(wf).Nt;
  else
    Nt = param.wfs(wf).Nt_raw;
  end
  Nc = size(param.load.imgs{img},1);
  data{img} = zeros(Nt,Nx,Nc,'single');
end

%% Load data
% ===================================================================
for board_idx = 1:length(boards)
  board = boards(board_idx);
  
  while rec = 1:total_rec;
    
    % Load record
    [state] = load_radar_data_fh(fns,offset,state);
    
    if ~param.proc.raw_data
      % Process all adc-wf pairs in this record
      for accum_idx = 1:length(state(board+1).wf)
        adc = state(board+1).adc(accum_idx);
        wf = state(board+1).wf(accum_idx);
        
        % Convert from quantization to voltage at the receiver input for the
        % maximum gain case:
        %  1. fast time gains less than the maximum for this record will be
        %     compensated for in the next step
        %  2. antenna effects not considered at this step
        state.data{adc,wf} = state.data{adc,wf} * wfs(wf).quantization_to_V * state.hdr{adc,wf}.quantization_to_V_adjustment;
        
        % Apply fast-time varying gain
        state.data{adc,wf} = state.data{adc,wf} .* interp1(wfs(wf).gain.Time, wfs(wf).gain.Gain, wfs(wf).time_raw(1:wfs(wf).Nt_raw));
      end
    end
    
    % Accumulate (presum)
    if adc<=size(tmp_data,1) && wf<=size(tmp_data,2) && ~isempty(tmp_data{adc,wf})
      if num_accum == 0
        state(board+1).data{accum_idx} = state.data{adc,wf};
      else
        state(board+1).data{accum_idx} = state(board+1).data{accum_idx} + tmp_data{adc,wf};
      end
    else
      if num_accum == 0
        state(board+1).data{accum_idx} = tmp_data{adc,wf};
      else
        state(board+1).data{accum_idx} = state(board+1).data{accum_idx} + tmp_data{adc,wf};
      end
    end
  end
  
  num_accum = num_accum + 1;
  if num_accum >= param.load.presums
    num_accum = 0;
    out_idx = out_idx + 1;
    for accum_idx = 1:length(state(board+1).wf)
      adc = state(board+1).adc(accum_idx);
      wf = state(board+1).wf(accum_idx);
      img = state(board+1).img(accum_idx);
      wf_adc_idx = state(board+1).wf_adc_idx(accum_idx);
      wf_adc_sum = state(board+1).wf_adc_sum(accum_idx);
      wf_adc_sum_done = state(board+1).wf_adc_sum_done(accum_idx);

      % Sum up wf-adc pairs until done
      state(board+1).data{accum_idx} = state(board+1).data{accum_idx} + 1i*sign(iq_mode)*state(board+1).data{accum_idx};
      if ~wf_adc_sum_done
        continue;
      end
 
  end
end

end







if param.load.file_version == 401
  load_radar_data_fh = @load_radar_data_401;
elseif param.load.file_version == 412
  load_radar_data_fh = @load_radar_data_412;
end
machineformat = 'ieee-be';
precision = 'int16';


fn = param.load.filenames{adc_idx}{fn_idx};

% Read up to block_size bytes of data
[fid,msg] = fopen(fn,'r',machineformat);



