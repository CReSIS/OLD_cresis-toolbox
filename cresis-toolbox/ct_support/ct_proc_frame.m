function proc_frm = ct_proc_frame(frm_type,frm_types_to_proc)
% proc_frm = ct_proc_frame(frm_type,frm_types_to_proc)
%
% Used by get_heights, csarp, etc. to determine if a frame should be
% processed or not.
%
% frm_type = The type of frame. Usually from frames file. frames.proc_mode(frame)
% frm_types_to_proc = 5 element cell array for:
%   R1: SAR/Doppler domain controls
%   R2: 0 = good frame/process, 1 = process but do not post, 2 = post only
%   R3: Reserved
%   R4: Reserved
%   U: User defined
%
% proc_frm = boolean indicating whether or not to process this frame
%
% Example:
% proc_frm = ct_proc_frame(10,{0,[0 10],0,0,0});

proc_frm = true;
for digit = 1:5
  frm_proc_mode_digit = mod(frm_type,10);
  frm_type = floor(frm_type/10);
  if ~any(frm_proc_mode_digit == frm_types_to_proc{digit})
    proc_frm = false;
    break;
  end
end

end
