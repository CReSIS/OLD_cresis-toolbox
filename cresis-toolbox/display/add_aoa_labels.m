function add_aoa_labels(ha,theta)
% function add_aoa_labels(ha,theta)
%
% Add angle of arrival x-axis labels for range-Doppler
%
% Author: Logan Smith

stheta = fftshift(theta);

add_theta = 0;
deg_labels = [];
while add_theta < max(theta)
  deg_labels(end+1) = add_theta;
  add_theta = add_theta + 10;
end
deg_labels(end+1) = max(theta);
deg_labels = [-fliplr(deg_labels(2:end)) deg_labels];
for tidx=1:length(deg_labels)
  if sign(deg_labels(tidx)) < 0
    xtick_locs(tidx) = find(stheta <= deg_labels(tidx),1,'last');
  else
    xtick_locs(tidx) = find(stheta >= deg_labels(tidx),1);
  end
end
set(ha,'XTick',xtick_locs)
set(ha,'XTickLabel',round(stheta(xtick_locs)))