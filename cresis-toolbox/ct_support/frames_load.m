function frames = frames_load(param)

if isempty(param)
  error('param is empty and should contain a radar parameter spreadsheet structure or a frames filename.');
end

if ischar(param)
  frames_fn = param
  param = [];
elseif isstruct(param)
  frames_fn = ct_filename_support(param,'','frames');
else
  error('param must be a string or struct.');
end

frames = load(frames_fn);
if ~isfield(frames,'frame_idxs')
  if isempty(param)
    if isfield(frames,'param')
      param = frames.param;
    else
      error('A param structure must be passed in instead of frames filename to use the old file format. After updating, the filename may be used with frames_load.');
    end
  end
  warning('Old frames file format. frames_update.m being run to update frames file.');
  frames_update(param,[]);
  frames = load(ct_filename_support(param,'','frames')); % Load "frames" variable
end
