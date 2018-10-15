function [hdr,data] = data_trim(hdr,data,time_trim)
% [hdr,data] = data_trim(hdr,data,time_trim)
%
% Author: John Paden

for img = 1:numel(data)
  Nt = length(hdr.time{img});
  dt = hdr.time{img}(2)-hdr.time{img}(1);
  fc = hdr.freq{img}(1);
  
  Nt = Nt - sum(time_trim);
  hdr.time{img} = hdr.time{img}(1+time_trim(1) : end-time_trim(2));
  T = Nt*dt;
  df = 1/T;
  hdr.freq{img} = fc + df * ifftshift(-floor(Nt/2) : floor((Nt-1)/2)).';
  
  data{img} = data{img}(time_trim(1) + (1:Nt),:,:);
end
