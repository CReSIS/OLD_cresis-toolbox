function [hdr,data] = data_load(param,records)
% [hdr,data] = data_load(param,records)

%% Preallocate data
% ===================================================================
[accum,state] = build_img_load_struct(param.load.imgs, param.load.adcs, ...
  struct('file_version',param.load.file_version));

total_rec = length(param.load.file_idx{1});
Nx = floor(total_rec/param.proc.presums);
if ~iscell(data)
  data = cell(size(param.load.imgs));
end
for img_idx = 1:length(param.load.imgs)
  wf = abs(param.load.imgs{img_idx}(1,1));
  if param.load.wf_adc_comb.en
    Nt = param.load.wf_adc_comb.Nt;
  else
    if param.proc.raw_data
      Nt = param.wfs(wf).Nt_raw;
    else
      Nt = param.wfs(wf).Nt;
    end
  end
  if param.proc.combine_rx
    Nc = 1;
  else
    Nc = size(param.load.imgs{img_idx},1);
  end
  if any([size(data{img_idx},1) size(data{img_idx},2) size(data{img_idx},3)] ~= [Nt Nx Nc])
    data{img_idx} = zeros(Nt,Nx,Nc,'single');
  elseif param.proc.combine_rx
    data{img_idx}(:) = 0;
  end
end

%% Load data
% ===================================================================
for board_idx = 1:length(boards)
  board = boards(board_idx);
  
  while rel_rec = 1:total_rec;
    
    % Load record
    [state] = load_radar_data_fh(fns,offset,state);
    
    % Process all adc-wf pairs in this record
    for accum_idx = 1:length(accum(board+1).wf)
      adc = accum(board+1).adc(accum_idx);
      wf = accum(board+1).wf(accum_idx);
      
      else
        % Convert from quantization to voltage at the receiver input for the
        % maximum gain case:
        %  1. fast time gains less than the maximum for this record will be
        %     undone in the next step
        %  2. antenna effects not considered at this step
        state.data{adc,wf} = state.data{adc,wf} * wfs(wf).quantization_to_V * state.hdr{adc,wf}.quantization_to_V_adjustment;
        
        % Apply fast-time varying gain
        state.data{adc,wf} = state.data{adc,wf} .* interp1(wfs(wf).gain.Time, wfs(wf).gain.Gain, wfs(wf).time_raw(1:wfs(wf).Nt_raw));
        old_tmp_data{adc,wf} = tmp_data{adc,wf};
    end
    
    % Accumulate (presum)
    if adc<=size(tmp_data,1) && wf<=size(tmp_data,2) && ~isempty(tmp_data{adc,wf})
      if num_accum == 0
        accum(board+1).data{accum_idx} = tmp_data{adc,wf};
      else
        accum(board+1).data{accum_idx} = accum(board+1).data{accum_idx} + tmp_data{adc,wf};
      end
    else
      if num_accum == 0
        accum(board+1).data{accum_idx} = tmp_data{adc,wf};
      else
        accum(board+1).data{accum_idx} = accum(board+1).data{accum_idx} + tmp_data{adc,wf};
      end
    end
  end
  
  num_accum = num_accum + 1;
  if num_accum >= param.proc.presums
    num_accum = 0;
    out_idx = out_idx + 1;
    for accum_idx = 1:length(accum(board+1).wf)
      adc = accum(board+1).adc(accum_idx);
      wf = accum(board+1).wf(accum_idx);
      img_idx = accum(board+1).img_idx(accum_idx);
      wf_adc_idx = accum(board+1).wf_adc_idx(accum_idx);
      iq_mode = accum(board+1).iq_mode(accum_idx);
      zero_pi_mode = accum(board+1).zero_pi_mode(accum_idx);
      
      % Combine I&Q channels if necessary
      if iq_mode == 1
        accum(board+1).data{accum_idx} = accum(board+1).data{accum_idx-1};
        continue;
      elseif abs(iq_mode) >= 2
        accum(board+1).data{accum_idx} = accum(board+1).data{accum_idx} + 1i*sign(iq_mode)*accum(board+1).data{accum_idx};
      end
      
      % Combine zero and pi channels if necessary
      if zero_pi_mode == 1
        accum(board+1).data{accum_idx} = accum(board+1).data{accum_idx-1};
        continue;
      elseif zero_pi_mode >= 2
        accum(board+1).data{accum_idx} = accum(board+1).data{accum_idx} + sign(zero_pi_mode)*accum(board+1).data{accum_idx};
      end
      
      % Store wf-adc pair in output matrix
      data{img_idx}(:,out_idx,wf_adc_idx) = accum(board+1).data{accum_idx} / param.proc.presums;
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



