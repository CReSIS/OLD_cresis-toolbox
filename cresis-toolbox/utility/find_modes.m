function modes = find_modes(pdf,num_modes)
% modes = find_modes(pdf,num_modes)
% 
% Finds the modes of a pdf or waveform up to a maximum of num_modes.
% The num_modes largest modes are returned.
%
% pdf = vector representing pdf or waveform
% num_modes = positive integer
% modes = the indices into the pdf vector
%
% Author: John Paden

modes = [];
upwards = true;
pdf = [pdf(:); pdf(1)];
for idx = 2:length(pdf)
  if pdf(idx) < pdf(idx-1) && upwards
    modes = [modes, idx-1];
    upwards = false;
  end
  if pdf(idx) > pdf(idx-1) && ~upwards
    upwards = true;
  end
end
[vals idxs] = sort(pdf(modes),'descend');
if length(idxs) < num_modes
  % Special case when not enough modes were found
  modes = modes([idxs ones(1,num_modes-length(idxs))]);
else
  modes = modes(idxs(1:num_modes));
end

return;
