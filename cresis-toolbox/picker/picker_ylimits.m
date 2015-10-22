function new_ylim = picker_ylimits(ylimit_str,cur_ylim)
% new_ylim = picker_ylimits(ylimit_str,cur_ylim)
%
% Support function for picker.m.  It is called from a couple
% places in picker_view and picker_pick to determine what
% the y-limits should be.
%
% Author: John Paden

try
  % Assumes the ylimit_str is a matlab expression that can be evaluated
  new_ylim = eval(sprintf('[%s]', ylimit_str));
catch
  fprintf('Y-Limits (us) field failed to evaluate\n');
  new_ylim = [-inf inf];
end

% Missing numbers are filled with inf
if isempty(new_ylim)
  new_ylim = [-inf inf];
elseif length(new_ylim) < 2
  new_ylim = [-inf new_ylim(1)];
end

% Inf is replaced with current y-limit bounds
if new_ylim(1) == -inf
  new_ylim(1) = cur_ylim(1);
end
if new_ylim(2) == inf
  new_ylim(2) = cur_ylim(2);
end

return;
