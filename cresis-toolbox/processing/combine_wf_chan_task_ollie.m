function combine_wf_chan_task_ollie(static_param_file_name,dynamic_param_file_name,frm,chunk_id)
% combine_wf_chan_task_ollie(static_param_file_name,dynamic_param_file_name,frm,chunk_id)
%
% This functions is compiled and used for manual job submission to Slurm on
% Ollie using batch_combine.sh.
%
if (ischar(frm))
  frm=str2num(frm);
end
if (ischar(chunk_id))
  chunk_id=str2num(chunk_id);
end

load(static_param_file_name,'static_param');
param=static_param;
load(dynamic_param_file_name,'dynamic_param');
param.load.frm = frm;
param.load.chunk_idx = chunk_id;
param.load.num_chunks = dynamic_param.frms.(['frm',num2str(frm)]).num_chunks;
if (frm>1)
  param.load.prev_frm_num_chunks = dynamic_param.frms.(['frm',num2str(frm-1)]).num_chunks;
else
  param.load.prev_frm_num_chunks = [];  
end

clear frm chunk_id;

if ~isfield(param.combine,'frm_types') || isempty(param.combine.frm_types)
  param.combine.frm_types = {-1,-1,-1,-1,-1};
end

% Remove frames that do not exist from param.cmd.frms list
load(ct_filename_support(param,'','frames')); % Load "frames" variable
if ~isfield(param.cmd,'frms') || isempty(param.cmd.frms)
  param.cmd.frms = 1:length(frames.frame_idxs);
end
[valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
if length(valid_frms) ~= length(param.cmd.frms)
  bad_mask = ones(size(param.cmd.frms));
  bad_mask(keep_idxs) = 0;
  warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
    param.cmd.frms(find(bad_mask,1)));
  param.cmd.frms = valid_frms;
end

if ~isfield(param.combine,'img_comb_layer_params')
  param.combine.img_comb_layer_params = [];
end

if ~isfield(param.combine,'trim_time')
  param.combine.trim_time = true;
end

if ~isfield(param.csarp,'pulse_comp') || isempty(param.csarp.pulse_comp)
  param.csarp.pulse_comp = 1;
end

if ~isfield(param.csarp,'presums') || isempty(param.csarp.presums)
  param.csarp.presums = 1;
end

if ~isfield(param.combine,'in_path') || isempty(param.combine.in_path)
  param.combine.in_path = 'out';
end

if ~isfield(param.combine,'array_path') || isempty(param.combine.array_path)
  param.combine.array_path = 'out';
end

if ~isfield(param.combine,'out_path') || isempty(param.combine.out_path)
  param.combine.out_path = param.combine.method;
end

if ~isfield(param.csarp,'out_path') || isempty(param.csarp.out_path)
  param.csarp.out_path = 'out';
end

if ~isfield(param.combine,'presums') || isempty(param.combine.presums)
  if ~isfield(param.csarp,'presums') || isempty(param.csarp.presums)
    param.combine.presums = 1;
  else
    param.combine.presums = param.csarp.presums;
  end
end

if ~isfield(param.combine,'sar_type') || isempty(param.combine.sar_type)
  if ~isfield(param.csarp,'sar_type') || isempty(param.csarp.sar_type)
    param.combine.sar_type = 'fk';
  else
    param.combine.sar_type = param.csarp.sar_type;
  end
end

if strcmpi(param.combine.sar_type,'f-k')
  error('Deprecated sar_type name. Change param.combine.sar_type from ''f-k'' to ''fk'' in  your parameters (or remove parameter since ''fk'' is the default mode).');
end

if ~isfield(param.combine,'chunk_len') || isempty(param.combine.chunk_len)
  if ~isfield(param.csarp,'chunk_len') || isempty(param.csarp.chunk_len)
    error('param.combine.chunk_len or param.csarp.chunk_len must be defined');
  else
    param.combine.chunk_len = param.csarp.chunk_len;
  end
end

% Handles multilooking syntax:
%  {{[1 1],[1 2],[1 3],[1 4],[1 5]},{[2 1],[2 2],[2 3],[2 4],[2 5]}}
%  If the image is a cell array it describes multilooking across apertures
if ~iscell(param.combine.imgs{1})
  % No special multilooking, reformat old syntax to new multilooking syntax
  for img = 1:length(param.combine.imgs)
    param.combine.imgs{img} = {param.combine.imgs{img}};
  end
end

for img = 1:length(param.combine.imgs)
  for ml_idx = 1:length(param.combine.imgs{img})
    % Imaginary image indices is for IQ combining during raw data load
    % which we do not need here.
    param.combine.imgs{img}{ml_idx} = abs(param.combine.imgs{img}{ml_idx});
  end
end

if ~isfield(param.combine,'out_path') || isempty(param.combine.out_path)
  param.combine.out_path = param.combine.method;
end

[success] = combine_wf_chan_task(param);