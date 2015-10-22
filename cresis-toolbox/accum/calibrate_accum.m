clear

base_path = '/mnt/tmp/20110329-accum/';

file_idxs = 178:185;

clear fdata weight;
for file_idx_idx = 1:length(file_idxs)
  file_idx = file_idxs(file_idx_idx);
  
  fn = get_filename(base_path,'accum_','',sprintf('%04d.dat',file_idx));
  fprintf('Loading file %d: %s\n', file_idx, fn);
  
  param.clk = 1e9/16;
  [hdr,data] = basic_load_accum(fn, param);

  % =======================================================================
  % Presum data
  presums = 10;
  num_rec = floor(size(data,2)/presums);
  for wf=1:size(data,3)
    data(:,1:num_rec,wf) = fir_dec(data(:,:,wf),presums);
  end
  data = data(:,1:num_rec,:);

  % =======================================================================
  % Convert from quantization to voltage @ ADC
  adc_bits = 14;
  Vpp_scale = 2;
  hw_presums = 8;
  bit_shifts = 2;
  data = (data-mean(data(:))) * Vpp_scale/2^adc_bits * 2^bit_shifts / hw_presums;

  prf_expected = 50e3/size(data,3)/hw_presums;

  Nt = size(data,1);
  fs = 1e9/8;
  df = fs/Nt;
  freq = (0:df:(Nt-1)*df).';

  % 30 MHz carrier frequency
  % Center frequencies: 600 to 900 MHz in 20 MHz steps
  fc = 875e6 - 20e6*(0:15);

  for wf = 1:size(data,3)
    fdata(:,file_idx_idx,wf) = mean(abs(fft(data(:,:,wf))).^2,2);
    weight(file_idx_idx,wf) = mean(abs(fdata(:,file_idx_idx,wf)).^2);
    fdata(1,file_idx_idx,wf) = -inf;
  end
end
weight_orig = weight;
fdata_orig = fdata;

fdata = fdata_orig;
weight = weight_orig;
weight(lp(weight)<-60) = 0;

clear final_psd;
for wf = 3:size(weight,2)
  for file_idx_idx = 1:size(weight,1)
    if weight(file_idx_idx,wf) > 0
      fdata(:,file_idx_idx,wf) = fdata(:,file_idx_idx,wf) / sqrt(weight(file_idx_idx,wf));
    else
      fdata(:,file_idx_idx,wf) = 0;
    end
  end
  final_psd(:,wf) = sum(fdata(:,:,wf),2) / sum(weight(:,wf) > 0);
end

clear final_weight
for wf = 3:size(weight,2)
  final_psd(:,wf) = final_psd(:,wf) ./ mean(final_psd(2:end,wf));
  final_weight(wf) = sum(weight(:,wf)) / sum(weight(:,wf) > 0);
end

for wf = 3:size(weight,2)
  wf
  for file_idx_idx = 1:size(weight,1)
    plot(fdata(:,file_idx_idx,wf));
    hold on;
  end
  plot(final_psd(:,wf),'r')
  hold off;
  pause;
end

save('/mnt/scratch2/mdce_tmp/accum_cal','final_psd','final_weight');

return;
