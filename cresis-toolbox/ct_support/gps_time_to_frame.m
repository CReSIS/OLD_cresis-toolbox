function frame_time = gps_time_to_frame(gps_time,frame_breaks)
% frame_time = gps_time_to_frame(gps_time,frame_breaks)
%
% Converts gps time into a psuedo frame time where frame time aligns with
% the framing of the segment. Frame time starts at 1 and ends at
% number_of_frames_in_segment+1.
%
% This function is useful for creating a slow time axis which is better
% connected to the frame structure (e.g. using GPS time in ANSI-C format is
% often not very useful for interpretation/analysis).
%
% gps_time: gps time vector from records file
% frame_breaks: one of two types:
%   'index' type: vector of frame indices from frames file
%     (i.e. frames.frame_idxs).
%   'gps_time' type: vector of gps times indicating the start gps time of
%     each frame.
% frame_break_type: string, default is 'index', but can be 'gps_time'
%  'index': Expects indices into gps_time to indicate frame starts
%  'gps_time': Expects gps times to indicate frame starts
%
% Example:
%   records_fn = '/cresis/snfs1/dataproducts/csarp_support/records/snow/2012_Greenland_P3/records_20120314_02.mat';
%   frames_fn = '/cresis/snfs1/dataproducts/csarp_support/frames/snow/2012_Greenland_P3/frames_20120314_02.mat';
%   records = load(records_fn);
%   load(frames_fn);
%   frame_time = gps_time_to_frame(records.gps_time,frames.frame_idxs);
%
% Author: John Paden

% If input frame_breaks is in gps time, convert to indices
if exist('frame_break_type','var') && strcmpi(frame_break_type,'gps_time')
  for idx = 1:length(frame_breaks)
    [~,frame_breaks(idx)] = min(gps_time - frame_breaks(idx));
  end
end
  
frm = 1;
if frm < length(frame_breaks)
  next_frm_gps_time = gps_time(frame_breaks(frm+1));
else
  next_frm_gps_time = gps_time(end) + median(diff(gps_time));
end

frame_time = zeros(size(gps_time));

for rec = 1:length(gps_time)
  % Determine which frame we are in if we are not at the last frame
  % Also determine the gps time of the start of the next frame
  if frm < length(frame_breaks)
    if rec >= frame_breaks(frm+1)
      frm = frm + 1;
      if frm < length(frame_breaks)
        next_frm_gps_time = gps_time(frame_breaks(frm+1));
      else
        next_frm_gps_time = gps_time(end) + median(diff(gps_time));
      end
    end
  end
  
  % Determine the fractional offset
  frm_fraction = (gps_time(rec) - gps_time(frame_breaks(frm))) ...
    ./ (next_frm_gps_time - gps_time(frame_breaks(frm)));
  
  frame_time(rec) = frm + frm_fraction;
end

return;
