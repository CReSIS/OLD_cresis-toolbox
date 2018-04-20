h_fig = figure(1); clf;
h_axes = axes;
hold(h_axes,'on');
xlim(h_axes,[-200 200]);
ylim(h_axes,[-100 100]);

y = [];
z = [];
for rline=1:size(array_param.fcs{1}{1}.pos,2)
  for idx=1:length(array_param.fcs{1})
    y(idx,rline) = array_param.fcs{1}{idx}.pos(2,rline);
    z(idx,rline) = array_param.fcs{1}{idx}.pos(3,rline);
  end
end
% imagesc(abs(y))
% colormap(jet);
% caxis([0 25]);

xlim(h_axes,[min(y(:)) max(y(:))]);
ylim(h_axes,[min(z(:)) max(z(:))]);


for idx=1:length(array_param.fcs{1})
  h_plot(idx) = plot(h_axes, NaN,NaN,'x');
end
h_title = title(h_axes, '');

for rline=1:10:size(array_param.fcs{1}{1}.pos,2)
  for idx=1:length(array_param.fcs{1})
    set(h_plot(idx),'XData',array_param.fcs{1}{idx}.pos(2,rline),'YData',array_param.fcs{1}{idx}.pos(3,rline));
  end
  set(h_title,'String',sprintf('%d',rline));
%   M(rline) = getframe;
  pause
end
% movie(M);
return


% For each trajectory find the sum of all the distances to other trajectories 
total_dists = zeros(length(array_param.fcs{1}),1);
for idx = 1:length(array_param.fcs{1})
  dist = bsxfun(@minus, y([1:idx-1 idx+1:end],:), y(idx,:)).^2 ...
    + bsxfun(@minus, z([1:idx-1 idx+1:end],:), z(idx,:)).^2;
  total_dists(idx) = sum(sum(dist));
end


