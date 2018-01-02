function frame_str_latex = ct_frame_latex(frame_str)
%
% Converts a frame string to a form that works well with tex and latex.
%
% frame_str: string containing the frame ID (e.g. '20141210_02_003')
%
% frame_str_latex: string with underscores 
%
% Example:
%  frame_str = '20141210_02_003';
%  title(ct_frame_latex(frame_str));
%  title(sprintf('%s: Echogram', ct_frame_latex(frame_str)));
%
% Author: John Paden

frame_str_latex = '';

for idx = 1:length(frame_str)
  if frame_str(idx) ~= '_'
    frame_str_latex(end+1) = frame_str(idx);
  else
    frame_str_latex(end+(1:2)) = ['\' frame_str(idx)];
  end
end
