function autogenerate_frames(param,param_override)
%
% Reads in a param file's generic column and creates frames for each
% segment selected.  Everything is automated.  This is useful for the
% FMCW radars where frames are always auto generated.
%
% param.records.frame_mode: decimal scalar interpretted as a bit mask
%  bit 1: autogenerate frames (i.e. call this function)
%  bit 2: special FMCW mode for breaking frames at settings changes
%  bit 3: overwrites frame without asking

% =====================================================================
% General Setup
% =====================================================================
fprintf('\n\n==============================================\n');

dbstack_info = dbstack;
if ~exist('param','var') || isempty(param) || length(dbstack_info) == 1
  % =====================================================================
  % Debug Setup
  % =====================================================================
  param = read_param_xls('/users/paden/scripts/branch/params-cr1/kuband2_param_2012_Greenland_P3.xls','20120423_01');
  
  clear('param_override');
  param_override.sched.type = 'no scheduler';
  
  % Input checking
  if ~exist('param','var')
    error('A struct array of parameters must be passed in\n');
  end
  global gRadar;
  if exist('param_override','var')
    param_override = merge_structs(gRadar,param_override);
  else
    param_override = gRadar;
  end
  
elseif ~isstruct(param)
  % Functional form
  param();
end
param = merge_structs(param, param_override);

fprintf('  autogenerate_frames %s\n', param.day_seg);
fprintf('==============================================\n\n');

if any(strcmpi(param.radar_name,{'icards','acords','mcords','mcrds','mcords2','mcords3','mcords4','mcords5'}))
  frame_length = 50000;
elseif any(strcmpi(param.radar_name,{'accum','accum2'}))
  frame_length = 20000;
elseif any(strcmpi(param.radar_name,{'kaband3','kuband','kuband2','kuband3','snow','snow2','snow3','snow5','snow8'}))
  frame_length = 5000;
else
  error('%s is not a supported radar', param.radar_name);
end
if isfield(param.records,'frame_length') && ~isempty(param.records.frame_length)
  frame_length = param.records.frame_length;
end
  
frame_mode = dec2bin(param.records.frame_mode,8);

fprintf('Autogenerating frames for %s (%s)\n', param.day_seg, datestr(now))

records_fn = ct_filename_support(param,'','records');
frames_fn = ct_filename_support(param,'','frames');

if exist(frames_fn,'file')
  fprintf('  Frame already exists %s\n', frames_fn);
  if frame_mode(end-2) == '0'
    user_input = input('  Type ''o'' to overwrite or return to skip: ', 's');
    if isempty(user_input) || ~strcmpi(user_input(1),'o')
      fprintf('    Skipping frame generation\n');
      return;
    end
  end
end

if ~exist(records_fn,'file')
  warning('  Skipping, records file does not exist %s', records_fn);
  return;
end

records = load(records_fn);
along_track = geodetic_to_along_track(records.lat,records.lon);

frame_breaks = 0:frame_length:along_track(end);
if along_track(end)-frame_breaks(end) < frame_length/2
  frame_breaks = frame_breaks(1:end-1);
end
frames = [];
frames.frame_idxs = zeros(size(frame_breaks));

frames.frame_idxs(1) = 1;
idx = 2;
rec = 2;
if any(strcmp(param.radar_name,{'snow2','kuband2'})) && frame_mode(end-1) == '1' ...
    && str2double(param.day_seg(1:8)) > 20120630
  % Only 2012 Antarctica DC8 and later snow2/kuband2 supports "settings" field, so the
  % "> 20120630" check supports this... NEEDS TO BE SET TO CORRECT DATE WHICH IS SOMETIME
  % IN THE MIDDLE OF THE CAMPAIGN...
  while idx <= length(frame_breaks)
    if along_track(rec) > frame_breaks(idx)
      frames.frame_idxs(idx) = rec;
      frames.nyquist_zone(idx) = double(mod(records.settings(rec),4));
      idx = idx + 1;
    elseif records.settings(rec) ~= records.settings(rec-1)
      if along_track(rec) - along_track(frames.frame_idxs(idx-1)) < frame_length/2 && idx > 2 && ~mod(floor(records.settings(rec-1)/2^2),2)
        % Remove last frame if it is short and not the first frame and
        % not a loopback segment
        idx = idx - 1;
      end
      new_frame_breaks = along_track(rec):frame_length:along_track(end);
      frame_breaks = [frame_breaks(1:idx-1) new_frame_breaks];
      if along_track(end)-frame_breaks(end) < frame_length/2 && idx ~= length(frame_breaks)
        frame_breaks = frame_breaks(1:end-1);
      end
      frames.frame_idxs(idx) = rec;
      frames.nyquist_zone(idx) = double(mod(records.settings(rec),4));
      idx = idx + 1;
    end
    rec = rec + 1;
  end
else
  while idx <= length(frame_breaks)
    if along_track(rec) > frame_breaks(idx)
      frames.frame_idxs(idx) = rec;
      idx = idx + 1;
    end
    rec = rec + 1;
  end
  frames.nyquist_zone = NaN*zeros(size(frames.frame_idxs));
end
frames.proc_mode = zeros(size(frames.frame_idxs));

fprintf('  Saving %s (%s)\n', frames_fn, datestr(now));
frames_fn_dir = fileparts(frames_fn);
if ~exist(frames_fn_dir,'dir')
  mkdir(frames_fn_dir);
end
save(frames_fn,'-v6','frames');

return;
