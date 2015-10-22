% script convert_records_to_new_format
%
% Convert from old records_YYYYMMDD_segS standard to 
% records_YYYYMMDD_SS standard
%
% Run with "run_convert_records_to_new_format_SEASON" scripts.

% =========================================================================
% Automated Work
% =========================================================================
fns = get_filenames(old_records_path,'records','seg','.mat');

for file_idx = 1:length(fns)
  [file_dir filename] = fileparts(fns{file_idx});
  filebase = filename(1:17);
  seg = str2double(filename(21:end));
  new_filename = sprintf('%s%02d.mat', filebase, seg);
  fprintf('Creating %s\n', new_filename);
  day_seg = sprintf('%s_%02d', filename(9:16), seg);
  param = [];
  for param_idx = 1:length(params)
    if strcmpi(params(param_idx).day_seg,day_seg)
      param = params(param_idx);
      break;
    end
  end
  if isempty(param)
    fprintf('  No parameter entry (skipping).\n');
    continue;
  else
    fprintf('  Found parameter entry.\n');
  end
  load(fns{file_idx});
  param_records = param.records;
  param_radar = param.radar;
  records = rmfield(records,'radar_config');
  records = rmfield(records,'param');
  records = rmfield(records,'notes');
  
  for adc_idx = 1:length(records.wfs)
    if isfield(records.wfs{adc_idx},'samp_delay')
      records.wfs{adc_idx}.t0 = records.wfs{adc_idx}.samp_delay - 11e-6;
      records.wfs{adc_idx} = rmfield(records.wfs{adc_idx},'samp_delay');
    end
  end
  
  for adc_idx = 1:length(records.filenames)
    for file_idx = 1:length(records.filenames{adc_idx})
      [tmp name ext] ...
        = fileparts(records.filenames{adc_idx}{file_idx});
      records.filenames{adc_idx}{file_idx} = [name ext];
    end
  end
  
  records_fn = fullfile(records_path,new_filename);
  save(records_fn,'records','param_records','param_radar');
  create_records_aux_files(records_fn);
end

