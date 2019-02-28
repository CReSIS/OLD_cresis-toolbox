function get_heights_task_ollie(static_param_file_name,dynamic_param_file_name,frm,break_id)
% get_heights_task_ollie(static_param_file_name,dynamic_param_file_name,frm,break_id)
%
% This functions is compiled and used for manual job submission to Slurm on
% Ollie using batch_qlook.sh.
%
if (ischar(frm))
  frm=str2num(frm);
end
if (ischar(break_id))
  break_id=str2num(break_id);
end

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

load(dynamic_param_file_name,'dynamic_param');
param.load.frm = dynamic_param.frms.(['frm',num2str(frm)]).frm_id;
param.load.recs = dynamic_param.frms.(['frm',num2str(frm)]).breaks.(['break',num2str(break_id)]).recs;
param.load.recs_keep = dynamic_param.frms.(['frm',num2str(frm)]).breaks.(['break',num2str(break_id)]).recs_keep;
param.load.imgs = param.get_heights.imgs;
clear frm break_id;

[success] = get_heights_task(param);