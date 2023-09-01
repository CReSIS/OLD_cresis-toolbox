function xo_report_figs(varargin)

global gRadar;

switch nargin
    case 1
        xo = varargin{1};
    case 2
        xo = varargin{1};
        xo_hdr = varargin{2};
end

if isempty(strfind(xo_hdr, '_'))
    save_figs = 0;
else
    save_figs = 1;
    reuse_loc = fullfile(gRadar.ct_tmp_path, 'radcal_mat', xo_hdr);
end

mean_hdr = sprintf('Mean xo angle %01.2f', mean(xo.cx_angle));
fig_hdr = [xo_hdr, ' ' mean_hdr];


%% geoplot
if 1

    figure;
    h_gx = geoaxes('Basemap','darkwater');
    geoplot(h_gx, xo.cx_lat, xo.cx_lon, 'x'); hold on;
    sgtitle([fig_hdr, ' : geoplot'], 'Interpreter', 'none');

    if save_figs
        set(findobj(gcf,'type','axes'),'FontWeight', 'Bold', 'FontSize',12);
        set(gcf, 'Position', get(0, 'Screensize'));
        fig_fn = fullfile(reuse_loc, sprintf('1report_geoplot_%s.fig', xo_hdr));
        fprintf('Saving %s\n', fig_fn);
        %   ct_saveas(gcf,fig_fn);
        fig_fn = fullfile(reuse_loc, sprintf('1report_geoplot_%s.png', xo_hdr));
        fprintf('Saving %s\n', fig_fn);
        ct_saveas(gcf,fig_fn);
    end
end

%% xo angles
if 1

    figure;

    min_val = floor(min(xo.cx_angle));
    max_val = ceil(max(xo.cx_angle));
    use_ylims = [ max([-3, min_val-1]), min([max_val+1,93]) ];

    subplot(121)
    plot(xo.cx_angle,'.');
    grid on; axis tight;
    %   ylim(use_ylims);
    %   yticks([0:10:90]);
    %   yticklabels([0:10:90]);
    xlabel('xo idx');
    ylabel('Angle, degrees');
    title('xo angle plot');

    subplot(122);
    histogram(xo.cx_angle, ...
        'DisplayStyle','bar', ... stairs
        'Orientation','horizontal', ...
        'FaceColor', '#D95319');
    grid on; axis tight;
    %   ylim(use_ylims);
    %   yticks([0:10:90]);
    %   yticklabels([0:10:90]);
    xlabel('Occurences');
    ylabel('Angle, degrees');
    title('xo angle distribution');

    sgtitle([fig_hdr, ' : Crossover Angles'], 'Interpreter', 'none');

    if save_figs
        set(findobj(gcf,'type','axes'),'FontWeight', 'Bold', 'FontSize',12);
        set(gcf, 'Position', get(0, 'Screensize'));
        fig_fn = fullfile(reuse_loc, sprintf('1report_angles_%s.fig', xo_hdr));
        fprintf('Saving %s\n', fig_fn);
        %   ct_saveas(gcf,fig_fn);
        fig_fn = fullfile(reuse_loc, sprintf('1report_angles_%s.png', xo_hdr));
        fprintf('Saving %s\n', fig_fn);
        ct_saveas(gcf,fig_fn);
    end

end

%% xo elev
if 1

    figure;
    subplot(121);
    grid on; hold on;
    plot(xo.pt1_elevation, 'r+');
    plot(xo.pt2_elevation, 'ro');
    xlabel('xo idx');
    ylabel('Elevation, m');
    legend('pt_1', 'pt_2');%, 'Interpreter', 'latex');

    offsets = xo.pt1_elevation-xo.pt2_elevation;

    subplot(122);
    grid on; hold on;
    bar(offsets)
    ylabel('abs elevation diff, m');
    xlabel('xo idx');
    if size(xo,1)<=20
        for idx_xo = 1:size(xo,1)
            xline(idx_xo, ':', ...
                sprintf('%s + %s', xo.pt1_frame{idx_xo}, xo.pt2_frame{idx_xo}), ...
                'LabelVerticalAlignment', 'top', 'LabelOrientation', 'aligned', ...
                'Interpreter', 'none');
        end
    end

    sgtitle([fig_hdr, ' : elevation'], 'Interpreter', 'none');

    if save_figs
        set(findobj(gcf,'type','axes'),'FontWeight', 'Bold', 'FontSize',12);
        set(gcf, 'Position', get(0, 'Screensize'));
        fig_fn = fullfile(reuse_loc, sprintf('1report_elev_%s.fig', xo_hdr));
        fprintf('Saving %s\n', fig_fn);
        %   ct_saveas(gcf,fig_fn);
        fig_fn = fullfile(reuse_loc, sprintf('1report_elev_%s.png', xo_hdr));
        fprintf('Saving %s\n', fig_fn);
        ct_saveas(gcf,fig_fn);
    end

end

%% 3D
if 1
    figure; hold on; grid on;
    plot(xo.cx_lon, xo.cx_lat, 'k.');
    plot3(xo.pt1_lon, xo.pt1_lat, xo.pt1_elevation, 'r+', 'LineWidth', 2);
    plot3(xo.pt2_lon, xo.pt2_lat, xo.pt2_elevation, 'ro', 'LineWidth', 2);

    plot3(xo.pt1_lon, xo.pt1_lat, xo.pt1_surface, 'g+', 'LineWidth', 2);
    plot3(xo.pt2_lon, xo.pt2_lat, xo.pt2_surface, 'go', 'LineWidth', 2);

    plot3(xo.pt1_lon, xo.pt1_lat, -xo.pt1_bottom, 'b+', 'LineWidth', 2);
    plot3(xo.pt2_lon, xo.pt2_lat, -xo.pt2_bottom, 'bo', 'LineWidth', 2);

    title([fig_hdr, ' : xo(lat,lon)[black.] RED[+o](elevation) GREEN[+o](surface) BLUE[+o](bottom)'] , 'Interpreter', 'none');

    if save_figs
        set(findobj(gcf,'type','axes'),'FontWeight', 'Bold', 'FontSize',12);
        set(gcf, 'Position', get(0, 'Screensize'));
        fig_fn = fullfile(reuse_loc, sprintf('1report_3D_%s.fig', xo_hdr));
        fprintf('Saving %s\n', fig_fn);
        %   ct_saveas(gcf,fig_fn);
        fig_fn = fullfile(reuse_loc, sprintf('1report_3D_%s.png', xo_hdr));
        fprintf('Saving %s\n', fig_fn);
        ct_saveas(gcf,fig_fn);
    end

end
