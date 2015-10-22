function doa_plot_estimates(param,results)

haxis_list = [];
for method = param.method.list
  for src_idx = 1:size(param.monte.DOA,2)
    h_fig = 100*src_idx + method;
    figure(h_fig); clf;
    haxis_list(end+1) = axes('Parent',h_fig);
    hold(haxis_list(end),'on');
    imagesc([],[],results.theta_est{method}(:,:,src_idx)*180/pi,'parent',haxis_list(end));
    h_color = colorbar;
    set(get(h_color,'YLabel'),'String','Theta ( \circ )');
    xlabel('Test')
    ylabel('Run')
    axis tight;
  end
end
linkaxes(haxis_list,'xy');

return
