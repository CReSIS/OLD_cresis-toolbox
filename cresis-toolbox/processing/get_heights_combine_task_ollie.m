function get_heights_combine_task_ollie(static_param_file_name)
% get_heights_combine_task_ollie(static_param_file_name)
%
% This functions is compiled and used for manual job submission to Slurm on
% Ollie using batch_qlook_2.sh.
%
load(static_param_file_name,'static_param');
param=static_param;

if isfield(param.get_heights,'qlook')
  warning('The get_heights.qlook field is deprecated. Please remove the qlook portion of each field (i.e. qlook.out_path should be just out_path).');
  qlook_fieldnames = fieldnames(param.get_heights.qlook);
  for name_idx = 1:length(qlook_fieldnames)
    param.get_heights.(qlook_fieldnames{name_idx}) = param.get_heights.qlook.(qlook_fieldnames{name_idx});
  end
  param.get_heights = rmfield(param.get_heights,'qlook');
end

success = get_heights_combine_task(param);