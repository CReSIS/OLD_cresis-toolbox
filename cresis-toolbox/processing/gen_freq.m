function [freq] = gen_freq(time,fc)

if ~exist('fc','var')
  fc = 0;
end

dt    = time(2) - time(1);
Nt    = length(time);
df    = 1/(Nt*dt);
freq  = fc + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';

return