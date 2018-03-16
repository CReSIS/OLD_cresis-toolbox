function plot_ascope(varargin)

start_idx = -160;
stop_idx = 40;
h_fig = figure(1);
clf(h_fig);
legend_str = {};
h_axes = axes('parent',h_fig);
for echo_idx = 1:2:length(varargin)
  for rline_idx = 1:length(varargin{echo_idx+1})
    rline = varargin{echo_idx+1}(rline_idx);
    [max_val,max_idx] = max(lp(varargin{echo_idx}(:,rline)));
    bins = max_idx + (start_idx:stop_idx);
    echo_vals = nan(size(bins));
    valid_bins = find(bins>=1 & bins<=size(varargin{echo_idx},1));
    echo_vals(valid_bins) = varargin{echo_idx}(bins(valid_bins),rline);
    xaxis = start_idx + (0:length(bins)-1);
    h_plot = plot(xaxis,lp(echo_vals)-max_val,'parent',h_axes);
    hold(h_axes,'on');
    grid(h_axes,'on');
    legend_str{end+1} = sprintf('%d:%d', (echo_idx-1)/2+1, rline);
  end
end
legend(h_axes,legend_str,'location','best');
xlim(xaxis([1 end]));
ylim([-80 0]);
hold(h_axes,'off');

end
