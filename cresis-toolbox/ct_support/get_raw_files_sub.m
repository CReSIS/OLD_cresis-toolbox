function load_info = get_raw_files_sub(param,records,record_idxs)



% adc_headers: the actual adc headers that were loaded
if ~isfield(param.records.file,'adc_headers') || isempty(param.records.file.adc_headers)
  param.records.file.adc_headers = param.records.file.adcs;
end

% boards_headers: the boards that the actual adc headers were loaded from
boards_headers = adc_to_board(param.radar_name,param.records.file.adcs);

for idx = 1:length(param.records.file.adc_headers)
  % adc: the specific ADC we would like to load_info
  adc = param.records.file.adc_headers(idx);
  % adc_idx: the records file index for this adc
  adc_idx = find(param.records.file.adcs == adc);
  if isempty(adc_idx)
    error('ADC %d not present in records file\n', adc);
  end
  
  % board: the board associated with the ADC we would like to load_info
  board = adc_to_board(param.radar_name,adc);
  % board_header: the board headers that we will use with this ADC
  board_header = adc_to_board(param.radar_name,param.records.file.adc_headers(adc_idx));
  % board_idx: the index into the records board list to use
  board_idx = find(board_header == boards_headers);
  
  % Just get the file-information for the records we need
  load_info.file_idx{idx} = relative_rec_num_to_file_idx_vector( ...
    record_idxs,records.relative_rec_num{board_idx});
  load_info.offset{idx} = records.offset(board_idx,:);
  file_idxs = unique(load_info.file_idx{idx});
  
  % Recognize if first record is really from previous file and it is a
  % valid record (i.e. offset does not equal -2^31)
  if sign(load_info.offset{idx}(1)) < 0 && load_info.offset{idx}(1) ~= -2^31
    file_idxs = [file_idxs(1)-1 file_idxs];
  end
  
  % Just copy the filenames we need
  load_info.filenames{idx}(file_idxs) = records.relative_filename{board_idx}(file_idxs);
  
  % Modify filename according to channel
  for file_idx = 1:length(load_info.filenames{idx})
    if ~isequal(param.records.file.adc_headers,param.records.file.adcs)
      load_info.filenames{idx}{file_idx}(9:10) = sprintf('%02d',board);
    end
  end
  
  filepath = get_segment_file_list(param,adc);
  
  % Convert relative file paths into absolute file paths if required,
  % also corrects filesep (\ and /)
  for file_idx = 1:length(load_info.filenames{idx})
    load_info.filenames{idx}{file_idx} ...
      = fullfile(filepath,load_info.filenames{idx}{file_idx});
  end
end
