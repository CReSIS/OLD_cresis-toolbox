hara;
% global gRadar;
flag_load_xo = 0;

%% initial setup
% seasons = eval(['{' sprintf('''201%d_Greenland_P3'' ', [4, 7:9]) '}']);
% seasons = eval(['{' sprintf('''201%d_Greenland_P3'' ', [0:9]) '}']);

%% specify table name to reuse mat file
% or leave it blank to load csv from OPS

%% 2009-2019

% seasons = {'2009_Greenland_P3', seasons{:}};

% xo_table_name = 'OPS_CReSIS_Crossovers_WKT_E7WRX9shd4';
% xo_angle_filter_str = 'le_0p02';
% xo_angle_filter_str = 'ge_89p9';

%% 2017

seasons = {'2017_Greenland_P3'};

xo_table_name = 'OPS_CReSIS_Crossovers_WKT_XYQQ1nEbE7';
% xo_angle_filter_str = 'le_0p5';
% xo_angle_filter_str = 'ge_85';
xo_angle_filter_str = 'ge_0';
% xo_angle_filter_str = 'le_0p239';


N_seasons = length(seasons);

%% Load xo

while flag_load_xo<3
    if ( ~exist('xo_table_name', 'var') || isempty(xo_table_name) )

        %% Setup OPS params

        ops_sys = 'rds';
        ops_param= [];
        ops_param.properties.location = 'arctic';
        if 1 % filter by season
            ops_param.properties.seasons = seasons;
        end

        % As it apears on ops.cresis.ku.edu for a WKT Projection WGS1984 (EPSG:4326)

        if 0
            % this is for the southweat coast (first trial)...
            ops_param.properties.bound = 'POLYGON((-41.50736205994954 63.094240466556094,-41.53020306022716 62.8157382852153,-41.00806137448067 62.80233791769747,-40.98097889949436 63.07202824018716,-41.50736205994954 63.094240466556094))';
            xo_table_name = 'OPS_CReSIS_Crossovers_WKT_uvixZqHCgm';

            % This is around the center....
            ops_param.properties.bound = 'POLYGON((-38.86555334533527 72.75036447802057,-37.47194556784986 72.46521344746948,-38.3466756631652 72.35393529569359,-38.965161734035966 72.49229739005605, -38.86555334533527 72.75036447802057))';
            xo_table_name = 'OPS_CReSIS_Crossovers_WKT_8Xvt4B5ejb'; % OPS_CReSIS_Crossovers_WKT_mlncYJyzUG

        elseif 0
            % non-coastal 2009-2019
            ops_param.properties.bound = 'POLYGON((-60.21185297752899 79.28631607719494,-38.990142955064044 80.77748577534513,-26.620239884631797 77.47703730992983,-31.313132039002255 72.9805972146493,-34.605965673152475 69.07115740728769,-42.40508975888023 65.73463547624439,-45.52641898905595 61.826899702579475,-47.29065655063459 65.37043635921154,-46.72929098881992 69.88001876512513,-51.721320562340466 73.25388326368206,-57.07865805987291 76.6739823649063,-61.56467548263641 77.77118856865765,-60.21185297752899 79.28631607719494))';
            xo_table_name = 'OPS_CReSIS_Crossovers_WKT_E7WRX9shd4';

        elseif 1
            % non-coastal 2017
            ops_param.properties.bound = 'POLYGON((-60.21185297752899 79.28631607719494,-38.990142955064044 80.77748577534513,-26.620239884631797 77.47703730992983,-31.313132039002255 72.9805972146493,-34.605965673152475 69.07115740728769,-42.40508975888023 65.73463547624439,-45.52641898905595 61.826899702579475,-47.29065655063459 65.37043635921154,-46.72929098881992 69.88001876512513,-51.721320562340466 73.25388326368206,-57.07865805987291 76.6739823649063,-61.56467548263641 77.77118856865765,-60.21185297752899 79.28631607719494))';
            xo_table_name = 'OPS_CReSIS_Crossovers_WKT_XYQQ1nEbE7';

        else
            %
            ops_param.properties.bound = 'POLYGON(())';

        end

        %% Query OPS

        if ( ~exist('xo_table_name', 'var') || isempty(xo_table_name) ) ...
                && ( ~exist('xo_table_tag', 'var') || isempty(xo_table_tag) )
            fprintf('OPS >> New Query >>\n');
            keyboard;
            [~, ops_csv] = opsGetCrossoversWithinPolygon( ops_sys, ops_param );
            ops_csv = ops_csv.properties;
            [~, xo_table_name] = fileparts(ops_csv);
        else
            fprintf('OPS >> REUSE CSV >>\n');
            ops_csv = ['data/csv/', xo_table_name, '.csv' ];
        end

        xo_table_tag = strrep(xo_table_name, 'OPS_CReSIS_Crossovers_WKT_', '');

        %% Load csv as table (setup non-OPS copy)

        try_ops = 1;

        ops_table_fn_tmp = fullfile(gRadar.ct_tmp_path, ...
            'radcal_mat', 'OPS_tables', [xo_table_name, '.mat']);

        % check if ct_tmp has a copy
        try
            if exist(ops_table_fn_tmp, 'file')
                fprintf('Loading XO table from ct_tmp: %s ', ops_table_fn_tmp);
                load(ops_table_fn_tmp);
                try_ops = 0;
                fprintf('Loaded\n');
            else
                fprintf('(no XO table in ct_tmp)\n');
            end
        catch
            fprintf('Failed\n');
        end

        % Now, try the OPS
        if try_ops
            ops_csv_fn = [gRadar.ops.url, ops_csv];
            fprintf('Loading CSV file from OPS: %s\n', ops_csv_fn);

            ops_weboptions = weboptions;
            ops_weboptions.CertificateFilename = '';
            ops_weboptions.ContentType = 'table';

            wait_or_quit = 1;
            iter = 100;
            each_wait = 5;
            fprintf('Waiting (%d secs) ', each_wait);
            while wait_or_quit
                try
                    xo_table_raw = webread(ops_csv_fn, ops_weboptions);
                    wait_or_quit = 0;
                    fprintf('Loaded');
                catch ME
                    fprintf('%d...', iter);
                    if strcmp(ME.identifier, 'MATLAB:webservices:ContentTypeMismatch')
                        pause(each_wait);
                        iter = iter-1;
                        wait_or_quit = wait_or_quit*iter;
                    else
                        fprintf('\n !!!! check ME');
                        wait_or_quit = 0;
                        %ops_weboptions.ContentType = 'text';
                        %csv_data = webread([gRadar.ops.url, ops_csv.properties], ops_weboptions);
                    end
                end
            end
            fprintf('\n');

            % save XO table in ct_tmp
            try
                fprintf('Saving XO table in ct_tmp: %s ', ops_table_fn_tmp);
                if ~exist(fileparts(ops_table_fn_tmp), 'dir')
                    mkdir(fileparts(ops_table_fn_tmp));
                end
                save(ops_table_fn_tmp, '-v7.3', 'xo_table_raw');
                fprintf('Saved\n');
            catch
                fprintf('Failed\n');
            end

            clear ops_sys ops_param ops_status ops_csv ops_csv_fn
            clear ops_weboptions wait_or_quit iter each_wait

        end


        %% Pre-Filter csv data

        xo = xo_table_raw;

        mask1 = zeros(size(xo_table_raw,1),1);
        mask2 = zeros(size(xo_table_raw,1),1);

        for idx_season = 1:N_seasons
            mask1 = [~cellfun(@isempty,regexp(xo.pt1_season, seasons(idx_season)))] + mask1;
            mask2 = [~cellfun(@isempty,regexp(xo.pt2_season, seasons(idx_season)))] + mask2;
        end

        good_idxs = find(mask1==mask2);
        xo = xo(good_idxs,:);

        if 0
            figure; hold on; grid on; % XMinorGrid ='on';
            plot(mask1,'x'); plot(mask2,'o');
            xticks(good_idxs);
            title('Matching xo points from given seasons');
        end
        clear mask1 mask2 idx_season good_idxs

        % xo_report_figs(xo,'Actual');

        % filter xo angles
        %%% good_idxs = find(xo.cx_angle <= 0.02);
        xo_angle_filter_operator = 'ge';
        xo_angle_filter_threshold = 0;

        % xo_angle_filter_operator = 'le';
        % xo_angle_filter_threshold = 0.239;

        eval( sprintf('good_idxs = find( %s(xo.cx_angle, %f) );', ...
            xo_angle_filter_operator,  xo_angle_filter_threshold) );

        xo = xo(good_idxs,:);

        % filter elev
        %     good_idxs = find(abs(xo.pt1_elevation-xo.pt2_elevation) <= 50);
        %     xo = xo(good_idxs,:);

        xo_angle_filter_str = [xo_angle_filter_operator, '_', ...
            char(strrep(string(xo_angle_filter_threshold), '.', 'p'))];
        xo_hdr_str = sprintf('%s_%s', xo_table_tag, xo_angle_filter_str);

        xo_report_figs(xo, xo_hdr_str );
        %     keyboard; % close all;

        %% Organize params{idx_season}(idx_day_seg)

        %N_seasons =

        % Do this later;

        %% Find rec and load data
        failed_idxs = [];

        N_xo = size(xo,1);
        for idx_xo = 1:N_xo
            ident = sprintf('%s_%s_xo_%04d', xo_table_tag, xo_angle_filter_str, idx_xo);
            fprintf('\n%s %d (of %d) >> qloox\n', ident, idx_xo, N_xo);
            try
                [qq(idx_xo)] = qloox(xo(idx_xo,:), ident);
            catch
                failed_idxs = [failed_idxs, idx_xo];
            end
        end

        failed_idxs2 = []
        if ~isempty(failed_idxs)
            N_xo_failed = length(failed_idxs);
            for idx_failed = 1:N_xo_failed
                idx_xo = failed_idxs(idx_failed);
                ident = sprintf('%s_%s_xo_%04d', xo_table_tag, xo_angle_filter_str, idx_xo);
                fprintf('\n%s %d (of %d) >> REDO qloox\n', ident, idx_failed, N_xo_failed);
                try
                    [qq(idx_xo)] = qloox(xo(idx_xo,:), ident);
                catch
                    failed_idxs2 = [failed_idxs2, idx_xo];
                end
            end
        end

        if N_xo == 0
            qq = [];
        end

        %% Save reuse_mat

        reuse_loc = fullfile(gRadar.ct_tmp_path, 'radcal_mat', ...
            [xo_table_tag, '_', xo_angle_filter_str]);
        reuse_fn = fullfile(reuse_loc, 'reuse.mat');
        fprintf('\n\nSaving REUSE MAT: %s >>', reuse_fn);
        if ~exist(reuse_loc, 'dir')
            mkdir(reuse_loc);
        end
        save(reuse_fn, '-v7.3', 'seasons', 'N_seasons', 'xo_table_raw', 'xo', 'N_xo', 'qq', 'xo_angle_filter_str', 'xo_hdr_str');
        fprintf('Saved\n');
        flag_load_xo = 10;

    else
        %% Load reuse_mat


        try
            xo_table_tag = strrep(xo_table_name, 'OPS_CReSIS_Crossovers_WKT_', '');
            reuse_loc = fullfile(gRadar.ct_tmp_path, 'radcal_mat', ...
                [xo_table_tag, '_', xo_angle_filter_str]);
            reuse_fn = fullfile(reuse_loc, 'reuse.mat');

            fprintf('Loading REUSE MAT: %s >>', reuse_fn);
            load(reuse_fn);
            fprintf('Loaded\n');
            flag_load_xo = 10;
        catch
            fprintf('\nERROR!!! NOT Loaded\n Trying again (%d/3)\n', flag_load_xo);
            flag_load_xo = flag_load_xo+1;
            if flag_load_xo == 3
                xo_table_name = []; % one final chance!
                flag_load_xo = -42;
                fprintf('Jump back to Querying OPS\n\n');
            end
        end
    end

end % while flag_load_xo

% return;

%% print gps times
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
fprintf('point\t\txo\t\t\t\t\t\trecord\t\t\t\t\tqlook\n');
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
for idx = 1:N_xo
    fprintf('xo idx %d\n', idx);
    fprintf('[1]gps:\t\t%f\t\t%f\t\t%f\n', ...
        xo.pt1_gps_time(idx), qq(idx).record{1}.gps_time, qq(idx).dd{1}.GPS_time);
    fprintf('[1]date:\t%s\t%s\t%s\n', ...
        gps2date( xo.pt1_gps_time(idx) ), ...
        gps2date( qq(idx).record{1}.gps_time), ...
        gps2date( qq(idx).dd{1}.GPS_time) );
    fprintf('[2]gps:\t\t%f\t\t%f\t\t%f\n', ...
        xo.pt2_gps_time(idx), qq(idx).record{2}.gps_time, qq(idx).dd{2}.GPS_time);
    fprintf('[2]date:\t%s\t%s\t%s\n', ...
        gps2date( xo.pt2_gps_time(idx) ), ...
        gps2date( qq(idx).record{2}.gps_time), ...
        gps2date( qq(idx).dd{2}.GPS_time) );
end
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
% pick_idxs = find(xo.cx_angle<0.239);

%% plots
figures_visible = 1;
xo_hdr_str = sprintf('%s_%s', xo_table_tag, xo_angle_filter_str);
fig_hdr = xo_hdr_str;

if 0 % redo qloox
    N_xo = size(xo,1);
    for idx_xo = 1:N_xo
        ident = sprintf('%s_%s_xo_%04d', xo_table_tag, xo_angle_filter_str, idx_xo);
        fprintf('\n%s %d (of %d) >> qloox\n', ident, idx_xo, N_xo);
        qloox(xo(idx_xo,:), ident);
    end
end

if 1
    %% elev axis
    h_fig_waveform_elev = figure('Name','waveform', 'visible', figures_visible);
    hold on;
    leg_str = [];

    try
        for idx_xo = 1:N_xo
            for idx_pt = 1:2
                plot(qq(idx_xo).dd{idx_pt}.elev_axis, 10*log10( qq(idx_xo).dd{idx_pt}.Data ) , '.-' );
                leg_str{(idx_xo-1)*2 + idx_pt} = eval(sprintf('xo.pt%d_frame{idx_xo}',idx_pt));
            end
        end
    catch
    end
    grid on;
    legend(leg_str, 'Interpreter','none');
    set(gca, 'Xdir', 'reverse');
    xlim([-4000 4000]);
    xlim([ 0 2500]);
    ylim([-210 -30]);
    xlabel('elev\_axis, meter');
    ylabel('Magnitude, dB');

    title(fig_hdr, 'Interpreter', 'None');
    fig_fn = fullfile(reuse_loc, sprintf('0summary_elev_%s.fig', fig_hdr));
    set(findobj(h_fig_waveform_elev,'type','axes'),'FontWeight', 'Bold', 'FontSize',12);
    set(h_fig_waveform_elev, 'Position', get(0, 'Screensize'));
    fprintf('Saving %s\n', fig_fn);
    ct_saveas(h_fig_waveform_elev,fig_fn);
    fig_fn = fullfile(reuse_loc, sprintf('0summary_elev_%s.png', fig_hdr));
    fprintf('Saving %s\n', fig_fn);
    ct_saveas(h_fig_waveform_elev,fig_fn);


end

if 1
    %% time axis
    h_fig_waveform_time = figure('Name','waveform2', 'visible', figures_visible);
    hold on;
    leg_str = [];

    try
        for idx_xo = 1:N_xo
            for idx_pt = 1:2
                plot(qq(idx_xo).dd{idx_pt}.time_axis/1e-6, 10*log10( qq(idx_xo).dd{idx_pt}.Data ) , '.-' );
                leg_str{(idx_xo-1)*2 + idx_pt} = eval(sprintf('xo.pt%d_frame{idx_xo}',idx_pt));
            end
        end
    catch
    end
    grid on;
    legend(leg_str, 'Interpreter','none');
    xlim([-15 55]);
    xlim([-8 6]);
    xlim([-4 4]);
    ylim([-140 -30]);
    xlabel('surface-zeroed Fast-Time, us');
    ylabel('Magnitude, dB');

    title(fig_hdr, 'Interpreter', 'None');
    fig_fn = fullfile(reuse_loc, sprintf('0summary_time_%s.fig', fig_hdr));
    set(findobj(h_fig_waveform_time,'type','axes'),'FontWeight', 'Bold', 'FontSize',12);
    set(h_fig_waveform_time, 'Position', get(0, 'Screensize'));
    fprintf('Saving %s\n', fig_fn);
    ct_saveas(h_fig_waveform_time,fig_fn);
    fig_fn = fullfile(reuse_loc, sprintf('0summary_time_%s.png', fig_hdr));
    fprintf('Saving %s\n', fig_fn);
    ct_saveas(h_fig_waveform_time,fig_fn);

end

if 1
    %% Pow deal
    P_act_vals = reshape([qq.P_act], [2, N_xo]).';
    P_est_vals = reshape([qq.P_est], [2, N_xo]).';

    h_fig_power = figure('Name','power', 'visible', figures_visible);

    subplot(211);
    b1 = bar(P_est_vals, 1, 'FaceAlpha', 0.21);
    hold on;
    b2 = bar(P_act_vals, 0.69, 'FaceAlpha', 1);
    b2(1).FaceColor = b1(1).FaceColor;
    b2(2).FaceColor = b1(2).FaceColor;

    set(gca, 'YDir', 'reverse');
    axis([0 N_xo+2 -90 -30]);
    xlabel('xo idx');
    ylabel('Power, dB');
    legend('P_est_pt1', 'p_est_pt2', 'P_act_pt1', 'P_act_pt2', ...
        'Interpreter', 'none', 'Location', 'SouthEast');
    title('Comparison of Power levels for two points at each xo');

    subplot(212);
    b = bar(diff(P_act_vals, 1, 2), 0.5, 'FaceColor', "#EDB120");
    axis([0 N_xo+2 -33 33]);
    xlabel('xo idx');
    ylabel('diff(Actual Power), dB');
    title('Difference in powers at each xo');

    sgtitle(fig_hdr, 'Interpreter', 'None');
    fig_fn = fullfile(reuse_loc, sprintf('0summary_power_%s.fig', fig_hdr));
    set(findobj(h_fig_power,'type','axes'),'FontWeight', 'Bold', 'FontSize',12);
    set(h_fig_power, 'Position', get(0, 'Screensize'));
    fprintf('Saving %s\n', fig_fn);
    ct_saveas(h_fig_power,fig_fn);
    fig_fn = fullfile(reuse_loc, sprintf('0summary_power_%s.png', fig_hdr));
    fprintf('Saving %s\n', fig_fn);
    ct_saveas(h_fig_power,fig_fn);

end
