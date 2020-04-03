function surfData_modify(param,fn,layers,varargin)
% tomo.surfData_modify(param,fn,layers,varargin)
%
% Function for modifying surfData files. Can be used to create new layers
% by specifying a layer that does not already exist.
%
% param: parameter spreadsheet structure
%  .radar_name: determines which radar to modify
%  .season_name: determines which season to modify
%  .day_seg: determines which segment to modify
%  .cmd.frms: determines which frames to modify
% fn: filename to ct_filename_out(), leave empty for default 'surfData'
% varargin: arbitrary list of name,value pairs.
%  third, fifth, etc. arguments are a string containing the name of the
%    property to modify
%  fourth, sixth, etc. arguments are the new value to use for the
%    corresponding named property
%  
% Author: John Paden

% Determine which frames to process
frames = frames_load(param);

if isempty(param.cmd.frms)
  param.cmd.frms = 1:length(frames.frame_idxs);
end
% Remove frames that do not exist from param.cmd.frms list
[valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
if length(valid_frms) ~= length(param.cmd.frms)
  bad_mask = ones(size(param.cmd.frms));
  bad_mask(keep_idxs) = 0;
  warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
    param.cmd.frms(find(bad_mask,1)));
  param.cmd.frms = valid_frms;
end

% out_path: output surfData directory
out_path = ct_filename_out(param, fn, 'CSARP_surfData');

for frm = param.cmd.frms
  fn = fullfile(out_path,sprintf('Data_%s_%03d.mat', param.day_seg, frm));
  
  % Load "surf" variable
  load(fn);
  for layer = layers(:).'
    for idx = 1:2:length(varargin)
      surf(layer).(varargin{idx}) = varargin{idx+1};
    end
  end
  save(fn,'-v7','surf');
end
