    %% Plot: 1. plot echogram
    h_fig_echo = figure(pass_idx); clf(h_fig_echo);
    set(h_fig_echo,'WindowStyle','docked')
    set(h_fig_echo,'NumberTitle','off')
    set(h_fig_echo,'Name',num2str(pass_idx))
    h_axes_echo(pass_idx) = axes('parent',h_fig_echo);
    imagesc([],pass(pass_idx).time*1e6,lp(pass(pass_idx).data),'parent', h_axes_echo(pass_idx)); %background image
    title_str = pass(pass_idx).param_pass.day_seg;
    title_str = regexprep(title_str,'_','\\_');
    title(h_axes_echo(pass_idx),title_str);
    colormap(h_axes_echo(pass_idx), 1-gray(256));
    hold(h_axes_echo(pass_idx),'on');
    if 1
      xlabel(h_axes_echo(pass_idx), 'Range line');
      ylabel(h_axes_echo(pass_idx), 'Two way travel time (\mus)');
    else
      [h_axes_echo_background,hp1,hp2] = plotyy(0:Nt-1,0:Nt-1,pass(pass_idx).time*1e6,pass(pass_idx).time*1e6,'parent',h_fig_echo);
      ylim(h_axes_echo_background(1),[0 Nt-1]);
      ylim(h_axes_echo_background(2),pass(pass_idx).time([1 end])*1e6);
      xlabel(h_axes_echo_background(1), 'Range line');
      ylabel(h_axes_echo_background(1), 'Range bin');
      ylabel(h_axes_echo_background(2), 'Time (\mus)');
      set(hp1,'XData',NaN,'YData',NaN);
      set(hp2,'XData',NaN,'YData',NaN);
      set(h_axes_echo_background(2),'YDir','reverse')
    end
%%    
h_fig = figure(999); % define figure handle
clf(h_fig); % clear figure 
h_axes = subplot(2,1,1,'parent',h_fig); % set handle for axes 1
h_axes(2) = subplot(2,1,2,'parent',h_fig);  % set handle for axes 2
h_axes = subplot(2,1,1,'parent',h_fig);
plot(1:10)
h_axes(2) = subplot(2,1,2,'parent',h_fig);
plot(10:20)
%h_axes = subplot(2,1,1,'parent',h_fig);
%subplot(2,1,1);
    
    %% Figure
    h_fig = fig
    
    figure(999);
    subplot(2,1,1)
    
    subplot(2,1,2)
    h1 = plot(AT_data.Btrack_End_Clip.P11/1e3, AT_data.elev_End_Clip.P2011);
    hold on
    h2 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip.P2014);
    h3 = plot(AT_data.Btrack_End_Clip.P18/1e3, AT_data.elev_End_Clip.P2018);
    title('Test Plot 1 - Alignment of Profiles is Correct');
    xlabel('Along Track (km)');
    ylabel('Elevation (m)');
    legend('original 2011', 'original 2014', 'original 2018', ...
     'Location', 'southeast');
    hold off
    
    line_locations = [500 1000 1500 2000 2500];
    line_locations_norm = interp1(xlims, cumsum(pos([1 3])), line_locations);
    
    for i=1:numel(line_locations_norm)
    annotation('line', ...
        [line_locations_norm(i) line_locations_norm(i)], ...
        [pos(2) pos(2)+pos(4)], ...
        'Color', 'r', ...
        'LineWidth', 1);
    end