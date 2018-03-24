% Get all the files for a whole season.
%
% Author: John Paden
%
% See also get_segment_file_list.m

params = read_param_xls('H:/scripts/ct_params/rds_param_2016_Greenland_TOdtu.xls');

for param_idx = 1:length(params)
  param = params(param_idx);
  [base_dir,adc_folder_name,fns,file_idxs] = get_segment_file_list(param,1,true);
  fprintf('%s\t%d\t%d\n', param.day_seg, file_idxs([1 end]))
end

return;
