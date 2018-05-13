if 1
%   records_fn = '/home/administrator/Scratch/csarp_support/records/rds/2018_Greenland_Polar6/records_20180510_01.mat';
%   surf_fn = '/home/administrator/Scratch/rds/2018_Greenland_Polar6/CSARP_noise/surf_20180510_01_img_01.mat';
  
  records_fn = '/home/administrator/Scratch/csarp_support/records/rds/2018_Greenland_Polar6/records_20180510_02.mat';
  surf_fn = '/home/administrator/Scratch/rds/2018_Greenland_Polar6/CSARP_noise/surf_20180510_02_img_01.mat';
% %   
  [surf_fn_dir,surf_fn_name,surf_fn_ext] = fileparts(surf_fn);
  archive2_fn = fullfile(surf_fn_dir,['archive2_', surf_fn_name, surf_fn_ext]);
  
  if ~exist(archive2_fn,'file')
    copyfile(surf_fn,archive2_fn)
  end
  
  bad_adcs = [1 2 3 4 5 8];
  
  if 1
    dd=load(surf_fn);
    max_vals = squeeze(max(dd.surf_vals));
  end
  
  Nx = size(dd.surf_vals,2);
  block = 1000;
  Nc = size(dd.surf_vals,3);
  ref_idx = 6;
  rlines = 1:block/2:Nx-block;
  
  lags = zeros(length(rlines),Nc);
  
  for adc_idx = 1:Nc
    for rline_idx = 1:length(rlines)
      rline = rlines(rline_idx);
      
      [acor,lag] = xcorr(max_vals(rline-1 + (1:block),adc_idx), max_vals(rline-1 + (1:block),ref_idx), 50);
      
      [~,max_idx] = max(acor);
      lags(rline_idx,adc_idx) = lag(max_idx);
      
    end
  end
  
end

[records_fn_dir,records_fn_name,records_fn_ext] = fileparts(records_fn);
archive2_fn = fullfile(records_fn_dir,['archive2_', records_fn_name, records_fn_ext]);

if ~exist(archive2_fn,'file')
  copyfile(records_fn,archive2_fn)
end

records = load(archive2_fn);

offsets = medfilt1(mean(lags(:,bad_adcs),2), 31);
offsets = medfilt1(round(offsets),101);

offsets = round(interp1(dd.gps_time(1,rlines), offsets, records.gps_time,'nearest','extrap'));
plot(offsets)


% Create file_num from records.relative_rec_num
for adc_idx = 1:length(records.relative_rec_num)
  for file_idx = 1:length(records.relative_rec_num{adc_idx})
    start = records.relative_rec_num{adc_idx}(file_idx);
    if file_idx == length(records.relative_rec_num{adc_idx})
      stop = length(records.gps_time);
    else
      stop = records.relative_rec_num{adc_idx}(file_idx+1)-1;
    end
    file_num{adc_idx}(start:stop) = file_idx;
  end
end

cur_offset = 0;
for rec = 1:size(records.offset,2)
  if offsets(rec) > cur_offset
    rec
    % Move records to the right
    shift = offsets(rec) - cur_offset;
    for adc = bad_adcs
      file_num{adc}(rec:end-shift) = file_num{adc}(rec+shift:end);
      file_num{adc}(end-shift+1:end) = file_num{adc}(end-shift);
    end
    records.offset(bad_adcs,rec:end-shift) = records.offset(bad_adcs,rec+shift:end);
    records.offset(bad_adcs,end-shift+1:end) = -2^31;
    cur_offset = offsets(rec);
  elseif offsets(rec) < cur_offset
    rec
    % Move records to the left
    shift = cur_offset - offsets(rec);
    for adc = bad_adcs
      file_num{adc}(rec+shift:end) = file_num{adc}(rec:end-shift);
      file_num{adc}(rec:rec+shift-1) = file_num{adc}(rec);
    end
    records.offset(bad_adcs,rec+shift:end) = records.offset(bad_adcs,rec:end-shift);
    records.offset(bad_adcs,rec:rec+shift-1) = -2^31;
    cur_offset = offsets(rec);
  end
end

% Update records.relative_rec_num from file_num
for adc = bad_adcs
  for file_idx = 1:length(records.relative_rec_num{adc})
    records.relative_rec_num{adc}(file_idx) = find(file_num{adc}==file_idx,1);
  end
end

save(records_fn,'-struct','records');
create_records_aux_files(records_fn);
