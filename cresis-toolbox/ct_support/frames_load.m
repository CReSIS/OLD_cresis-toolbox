function frames = frames_load(param)

frames_fn = ct_filename_support(param,'','frames');
frames = load(frames_fn);
if ~isfield(frames,'frame_idxs')
  warning('Old frames file format. frames_update.m being run to update frames file.');
  frames_update(param,[]);
  frames = load(ct_filename_support(param,'','frames')); % Load "frames" variable
end
