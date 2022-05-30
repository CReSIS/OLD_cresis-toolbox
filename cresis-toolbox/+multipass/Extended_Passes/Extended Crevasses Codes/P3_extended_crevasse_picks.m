%% clip interpolated datasets to grounding line location
for i = 1:numel(XY_data)
    AT_data(i).Interp_bed = XY_data(i).Interp_Bed_Corrected(419:3083);
    AT_data(i).Interp_AT_dist = XY_data(i).Interp_AT_Dist(419:3083);
    AT_data(i).Interp_AT_dist_0 = AT_data(i).Interp_AT_dist - AT_data(i).Interp_AT_dist(1);
    AT_data(i).InterP_vel_mag = XY_data(i).Interp_Vel_Mag(419:3083);
    AT_data(i).Interp_AT_vel = AT_data(i).Interp_AT_dist + AT_data(i).InterP_vel_mag;
    AT_data(i).Interp_AT_vel_0 = AT_data(i).Interp_AT_dist - AT_data(i).Interp_AT_dist(1);
end

%% Basal Crevasse Apex Picking section:
% Crevasse Apex calculation and binary file formation
for i = 1:numel(XY_data)
    AT_data(i).apex_bin = islocalmax(AT_data(i).Interp_bed);
    AT_data(i).apex_pic = AT_data(i).Interp_bed.*AT_data(i).apex_bin;
    AT_data(i).apex_pic(AT_data(i).apex_pic == 0) = NaN;
    AT_data(i).index_apex = find(~isnan(AT_data(i).apex_pic));
    AT_data(i).AT_dist_apex_vel = AT_data(i).Interp_AT_vel(AT_data(i).index_apex);
    AT_data(i).AT_dist_apex = AT_data(i).Interp_AT_dist(AT_data(i).index_apex);
    AT_data(i).crevasse_apex = AT_data(i).Interp_bed(AT_data(i).index_apex);
end
    
%% Basal Crevasse Wall Base Picking:
% Crevasse wall base calculation and binary file formation
for i = 1:numel(XY_data)
    AT_data(i).base_bin = islocalmin(AT_data(i).Interp_bed);
    AT_data(i).base_pic = AT_data(i).Interp_bed.*AT_data(i).base_bin;
    AT_data(i).base_pic(AT_data(i).base_pic == 0) = NaN;
    AT_data(i).index_base = find(~isnan(AT_data(i).base_pic));
    AT_data(i).AT_dist_base_vel = AT_data(i).Interp_AT_vel(AT_data(i).index_base);
    AT_data(i).AT_dist_base = AT_data(i).Interp_AT_dist(AT_data(i).index_base);
    AT_data(i).crevasse_base =  AT_data(i).Interp_bed(AT_data(i).index_base);
end

%% Basal Crevasse Apex Picks
% short side denoted in binary array direction_bin, r=1, l=0
AT_data(1).apexes_idx = [4, 12, 23, 40, 44, 46, 62, 64, 82, 93, 110, 122, 137];
%AT_data(1).apex_direction = [r, l, r, r, r, r, l, r, r, l];
% Get index of left and right short sides
AT_data(1).apex_direction_bin_r = [1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1];
AT_data(1).apex_direction_bin_l = [0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0];
AT_data(1).apex_data = AT_data(1).crevasse_apex(AT_data(1).apexes_idx);
AT_data(1).base_data = AT_data(1).crevasse_base(AT_data(1).apexes_idx);
AT_data(1).apex_AT = AT_data(1).AT_dist_apex(AT_data(1).apexes_idx);
AT_data(1).left = [2, 12, 21, 39, 43, 45, 61, 62, 80, 91, 109, 118, 133]; 
AT_data(1).left_base = AT_data(1).AT_dist_base(AT_data(1).left);
AT_data(1).right = [5, 13, 25, 42, 44, 46, 62, 65, 82, 93, 111, 126, 138];
AT_data(1).right_base = AT_data(1).AT_dist_base(AT_data(1).right);
AT_data(1).year_array = 2010*ones(length(AT_data(1).apexes_idx),1);

AT_data(2).apexes_idx = [8, 20, 33, 62, 69, 73, 97, 104, 130, 145, 169, 177, 201];
%AT_data(2).apex_direction = [r, r, l, r, l, l, r, r, l, l];
AT_data(2).apex_direction_bin_r = [1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0];
AT_data(2).apex_direction_bin_l = [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1];
AT_data(2).apex_data = AT_data(2).crevasse_apex(AT_data(2).apexes_idx);
AT_data(2).apex_AT = AT_data(2).AT_dist_apex(AT_data(2).apexes_idx);
AT_data(2).left = [4, 19, 30, 58, 68, 71, 94, 99, 127, 141, 168, 176, 193]; 
AT_data(2).left_base = AT_data(2).AT_dist_base(AT_data(2).left);
AT_data(2).right = [9, 21, 36, 63, 70, 74, 99, 107, 131, 148, 173, 181, 202];
AT_data(2).right_base = AT_data(2).AT_dist_base(AT_data(2).right);
AT_data(2).year_array = 2011*ones(length(AT_data(2).apexes_idx),1);

AT_data(3).apexes_idx = [9, 22, 33, 63, 72, 76, 92, 99, 123, 143, 167, 176, 205];
%AT_data(3).apex_direction = [r, r, l, r, l, l, r, r, l, l];
AT_data(3).apex_direction_bin_r = [1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0];
AT_data(3).apex_direction_bin_l = [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1];
AT_data(3).apex_data = AT_data(2).crevasse_apex(AT_data(2).apexes_idx);
AT_data(3).apex_AT = AT_data(2).AT_dist_apex(AT_data(2).apexes_idx);
AT_data(3).left = [9, 22, 30, 62, 69, 74, 92, 96, 120, 138, 167, 175, 202]; 
AT_data(3).left_base = AT_data(2).AT_dist_base(AT_data(2).left);
AT_data(3).right = [10, 22, 34, 65, 74, 78, 95, 102, 124, 144, 173, 187, 207];
AT_data(3).right_base = AT_data(2).AT_dist_base(AT_data(2).right);
AT_data(3).year_array = 2013*ones(length(AT_data(3).apexes_idx),1);

AT_data(4).apexes_idx = [20, 27, 37, 76, 83, 92, 127, 139, 173, 196, 223, 236, 277];
%AT_data(4).apex_direction = [r, r, l, r, l, l, r, r, l, l];
AT_data(4).apex_direction_bin_r = [1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1];
AT_data(4).apex_direction_bin_l = [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0];
AT_data(4).apex_data = AT_data(4).crevasse_apex(AT_data(4).apexes_idx);
AT_data(4).apex_AT = AT_data(4).AT_dist_apex(AT_data(4).apexes_idx);
AT_data(4).left = [16, 26, 36, 71, 79, 90, 123, 129, 166, 190, 222, 235, 266]; 
AT_data(4).left_base = AT_data(4).AT_dist_base(AT_data(4).left);
AT_data(4).right = [22, 27, 40, 79, 90, 94, 130, 142, 174, 196, 232, 239, 278];
AT_data(4).right_base = AT_data(4).AT_dist_base(AT_data(4).right);
AT_data(4).year_array = 2014*ones(length(AT_data(4).apexes_idx),1);

AT_data(5).apexes_idx = [7, 20, 30, 60, 66, 66, 95, 100, 121, 137, 156, 168, 181];
%AT_data(5).apex_direction = [r, r, l, r, l, l, r, r, l, l];
AT_data(5).apex_direction_bin_r = [1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1];
AT_data(5).apex_direction_bin_l = [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0];
AT_data(5).apex_data = AT_data(5).crevasse_apex(AT_data(5).apexes_idx);
AT_data(5).apex_AT = AT_data(5).AT_dist_apex(AT_data(5).apexes_idx);
AT_data(5).left = [7, 18, 27, 54, 62, 62, 91, 97, 115, 129, 156, 165, 179]; 
AT_data(5).left_base = AT_data(5).AT_dist_base(AT_data(5).left);
AT_data(5).right = [16, 21, 31, 62, 72, 72, 97, 101, 122, 132, 160, 176, 182];
AT_data(5).right_base = AT_data(5).AT_dist_base(AT_data(5).right);
AT_data(5).year_array = 2017*ones(length(AT_data(5).apexes_idx),1);

for i = 1:5
    AT_data(i).apex_data = AT_data(i).crevasse_apex(AT_data(i).apexes_idx);
    %AT_data(i).base_data = AT_data(i).crevasse_base(AT_data(i).apexes_idx);
    AT_data(i).base_width = AT_data(i).right_base - AT_data(i).left_base;
    AT_data(i).short_idx_r = find(AT_data(i).apex_direction_bin_r == 1);
    AT_data(i).short_idx_l = find(AT_data(i).apex_direction_bin_l == 1);
    AT_data(i).short_idx_r_pos = AT_data(i).right(AT_data(i).short_idx_r);
    AT_data(i).short_idx_l_pos = AT_data(i).left(AT_data(i).short_idx_l);
    AT_data(i).short_base(AT_data(i).short_idx_r) = (AT_data(i).crevasse_base(AT_data(i).short_idx_r_pos)); 
    AT_data(i).short_base(AT_data(i).short_idx_l) = (AT_data(i).crevasse_base(AT_data(i).short_idx_l_pos));
    AT_data(i).long_idx_r = find(AT_data(i).apex_direction_bin_r == 0);
    AT_data(i).long_idx_l = find(AT_data(i).apex_direction_bin_l == 0);
    AT_data(i).long_base = [AT_data(i).base_data(AT_data(i).long_idx_r) AT_data(i).base_data(AT_data(i).long_idx_l)];
    AT_data(i).DH_long = AT_data(i).apex_data - AT_data(i).long_base;
    AT_data(i).DH_short = AT_data(i).apex_data - AT_data(i).short_base;
    AT_data(i).AT_apex_0km = AT_data(i).AT_dist_apex_vel/1e3 - AT_data(i).Interp_AT_vel(1)/1e3;
    AT_data(i).AT_base_0km = AT_data(i).AT_dist_base_vel/1e3 - AT_data(i).Interp_AT_vel(1)/1e3;
    AT_data(i).AT_interp_0km = AT_data(i).Interp_AT_vel/1e3 - AT_data(i).Interp_AT_vel(1)/1e3;
end

% Define time period off summations for annual correction
time_offset.Period_1 = time_offset.P10 + time_offset.P11;
time_offset.Period_2 = time_offset.P11_2 + time_offset.P12 + time_offset.P13;
time_offset.Period_3 = time_offset.P13_2 + time_offset.P14;
time_offset.Period_4 = time_offset.P14_2 + time_offset.P15 + time_offset.P16 + time_offset.P17;

% loop for taking manual data and splitting into groupings per crevasse
for i = 1:length(AT_data(5).apexes_idx)
    plot_data(i).short_H = [AT_data(1).DH_short(i), AT_data(2).DH_short(i), AT_data(3).DH_short(i),...
        AT_data(4).DH_short(i),AT_data(5).DH_short(i)];
    plot_data(i).year_array = [AT_data(1).year_array(i), AT_data(2).year_array(i),...
        AT_data(3).year_array(i), AT_data(4).year_array(i), AT_data(5).year_array(i)];
    plot_data(i).apex_elev = [AT_data(1).apex_data(i), AT_data(3).apex_data(i), ...
        AT_data(3).apex_data(i), AT_data(4).apex_data(i), AT_data(5).apex_data(i)];
    plot_data(i).apex_GL_dist = [AT_data(1).apex_AT(i)/1e3, AT_data(2).apex_AT(i)/1e3,...
        AT_data(3).apex_AT(i)/1e3, AT_data(4).apex_AT(i)/1e3, AT_data(5).apex_AT(i)/1e3];
    plot_data(i).DH_short = [AT_data(1).DH_short(i), AT_data(2).DH_short(i), ...
        AT_data(3).DH_short(i), AT_data(4).DH_short(i), AT_data(5).DH_short(i)];
    plot_data(i).DH_long = [AT_data(1).DH_long(i), AT_data(2).DH_long(i), ...
        AT_data(3).DH_long(i), AT_data(4).DH_long(i), AT_data(5).DH_long(i)];
    plot_data(i).widths = [AT_data(1).base_width(i), AT_data(2).base_width(i), ...
        AT_data(3).base_width(i), AT_data(4).base_width(i), AT_data(5).base_width(i)];
    plot_data(i).DW_annual = [((AT_data(2).base_width(i) - AT_data(1).base_width(i))/(time_offset.Period_1)), ...
        ((AT_data(3).base_width(i) - AT_data(2).base_width(i))/(time_offset.Period_2)),...
        ((AT_data(4).base_width(i) - AT_data(3).base_width(i))/(time_offset.Period_3)),...
        ((AT_data(5).base_width(i) - AT_data(4).base_width(i))/(time_offset.Period_4))];
    plot_data(i).DH_short_annual = [((AT_data(2).DH_short(i) - AT_data(1).DH_short(i))/(time_offset.Period_1)), ...
        ((AT_data(3).DH_short(i) - AT_data(2).DH_short(i))/(time_offset.Period_2)),...
        ((AT_data(4).DH_short(i) - AT_data(3).DH_short(i))/(time_offset.Period_3)),...
        ((AT_data(5).DH_short(i) - AT_data(4).DH_short(i))/(time_offset.Period_4))];
    plot_data(i).DH_long_annual = [((AT_data(2).DH_long(i) - AT_data(1).DH_long(i))/(time_offset.Period_1)), ...
        ((AT_data(3).DH_long(i) - AT_data(2).DH_long(i))/(time_offset.Period_2)),...
        ((AT_data(4).DH_long(i) - AT_data(3).DH_long(i))/(time_offset.Period_3)),...
        ((AT_data(5).DH_long(i) - AT_data(4).DH_long(i))/(time_offset.Period_4))];
    plot_data(i).DX = [(AT_data(2).apex_AT(i)/1e3 - AT_data(1).apex_AT(i)/1e3), (AT_data(3).apex_AT(i)/1e3 - AT_data(2).apex_AT(i)/1e3),...
        (AT_data(4).apex_AT(i)/1e3 - AT_data(3).apex_AT(i)/1e3), (AT_data(5).apex_AT(i)/1e3 - AT_data(4).apex_AT(i)/1e3)];
    plot_data(i).DX_yr = [((AT_data(2).apex_AT(i)/1e3 - AT_data(1).apex_AT(i)/1e3)/time_offset.Period_1),...
        ((AT_data(3).apex_AT(i)/1e3 - AT_data(2).apex_AT(i)/1e3)/time_offset.Period_2),...
        ((AT_data(4).apex_AT(i)/1e3 - AT_data(3).apex_AT(i)/1e3)/time_offset.Period_3),...
        ((AT_data(5).apex_AT(i)/1e3 - AT_data(4).apex_AT(i)/1e3)/time_offset.Period_4)];
    plot_data(i).DX_midpoint = [((AT_data(2).apex_AT(i)/1e3 + AT_data(1).apex_AT(i)/1e3)/2), ...
        ((AT_data(3).apex_AT(i)/1e3 + AT_data(2).apex_AT(i)/1e3)/2),...
        ((AT_data(4).apex_AT(i)/1e3 + AT_data(3).apex_AT(i)/1e3)/2),...
        ((AT_data(5).apex_AT(i)/1e3 + AT_data(4).apex_AT(i)/1e3)/2)];
    plot_data(i).ratio_H_W = [AT_data(1).DH_short(i)/AT_data(1).base_width(i), AT_data(2).DH_short(i)/AT_data(2).base_width(i),...
        AT_data(3).DH_short(i)/AT_data(3).base_width(i), AT_data(4).DH_short(i)/AT_data(4).base_width(i),...
        AT_data(5).DH_short(i)/AT_data(5).base_width(i)];
    plot_data(i).ratio_W_H = [AT_data(1).base_width(i)/AT_data(1).DH_short(i), AT_data(2).base_width(i)/AT_data(2).DH_short(i),...
        AT_data(3).base_width(i)/AT_data(3).DH_short(i), AT_data(4).base_width(i)/AT_data(4).DH_short(i),...
        AT_data(5).base_width(i)/AT_data(5).DH_short(i)];
    plot_data(i).mean_H_W_ratio = mean(plot_data(i).ratio_H_W);
    plot_data(i).std_H_W_ratio = std(plot_data(i).ratio_H_W);
end

% mean
plot_data(i).mean_HW_ratio_total = mean([plot_data(1).mean_H_W_ratio, plot_data(2).mean_H_W_ratio, plot_data(3).mean_H_W_ratio,...
    plot_data(4).mean_H_W_ratio, plot_data(5).mean_H_W_ratio, plot_data(6).mean_H_W_ratio, plot_data(7).mean_H_W_ratio, ...
    plot_data(8).mean_H_W_ratio, plot_data(9).mean_H_W_ratio, plot_data(10).mean_H_W_ratio, plot_data(11).mean_H_W_ratio, ...
    plot_data(12).mean_H_W_ratio, plot_data(13).mean_H_W_ratio]);
%std
plot_data(i).std_HW_ratio_total = std([plot_data(1).std_H_W_ratio, plot_data(2).std_H_W_ratio, plot_data(3).std_H_W_ratio,...
    plot_data(4).std_H_W_ratio, plot_data(5).std_H_W_ratio, plot_data(6).std_H_W_ratio, plot_data(7).std_H_W_ratio, ...
    plot_data(8).std_H_W_ratio, plot_data(9).std_H_W_ratio, plot_data(10).std_H_W_ratio, plot_data(11).std_H_W_ratio, ...
    plot_data(12).std_H_W_ratio, plot_data(13).std_H_W_ratio]);

for i = 1:length(AT_data(5).apexes_idx)
    plot_data(i).normalized_HW_ratio = [((plot_data(i).ratio_H_W(1) - plot_data(1).mean_HW_ratio_total(1))/plot_data(1).std_HW_ratio_total(1)),... 
    ((plot_data(i).ratio_H_W(2) - plot_data(1).mean_HW_ratio_total(1))/plot_data(1).std_HW_ratio_total(1)), ...
    ((plot_data(i).ratio_H_W(3) - plot_data(1).mean_HW_ratio_total(1))/plot_data(1).std_HW_ratio_total(1)), ...
    ((plot_data(i).ratio_H_W(4) - plot_data(1).mean_HW_ratio_total(1))/plot_data(1).std_HW_ratio_total(1)), ...
    ((plot_data(i).ratio_H_W(5) - plot_data(1).mean_HW_ratio_total(1))/plot_data(1).std_HW_ratio_total(1))];
end

%% Figures for Poster
   % Basal Epoch Change with respect to mid point 
colors = colormap(flipud(turbo(length(plot_data))));
figure(8)
hold on 
for k = 1:numel(plot_data)
    legend_text =['Crevasse ', num2str(k)];
    size = 90;
    plot(plot_data(k).year_array, plot_data(k).ratio_W_H,'-o', 'MarkerFaceColor',colors(k,:),...
        'Color',colors(k,:), 'DisplayName', legend_text);
    h = plot(plot_data(k).year_array, plot_data(k).ratio_W_H, 'o', 'Color','k');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h1 = scatter(plot_data(k).year_array(1), plot_data(k).ratio_W_H(1), size,'o', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h2 = scatter(plot_data(k).year_array(2), plot_data(k).ratio_W_H(2), size,'v', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h3 = scatter(plot_data(k).year_array(3), plot_data(k).ratio_W_H(3), size,'s', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h3.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h4 = scatter(plot_data(k).year_array(4), plot_data(k).ratio_W_H(4), size,'p', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h5 = scatter(plot_data(k).year_array(5), plot_data(k).ratio_W_H(5), size,'h', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h5.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
   title('P3 Change in Crevasse Width to Height Ratio Per Epoch');
   xlabel('Year');
   ylabel('Crevasse H/W');
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;

figure(9)
hold on 
for k = 1:numel(plot_data)
    legend_text =['Crevasse ', num2str(k)];
    size = 90;
    plot(plot_data(k).year_array, plot_data(k).ratio_H_W,'-o', 'MarkerFaceColor',colors(k,:),...
        'Color',colors(k,:), 'DisplayName', legend_text);
    h = plot(plot_data(k).year_array, plot_data(k).ratio_H_W, 'o', 'Color','k');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h1 = scatter(plot_data(k).year_array(1), plot_data(k).ratio_H_W(1), size,'o', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h2 = scatter(plot_data(k).year_array(2), plot_data(k).ratio_H_W(2), size,'v', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h3 = scatter(plot_data(k).year_array(3), plot_data(k).ratio_H_W(3), size,'s', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h3.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h4 = scatter(plot_data(k).year_array(4), plot_data(k).ratio_H_W(4), size,'p', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h5 = scatter(plot_data(k).year_array(5), plot_data(k).ratio_H_W(5), size,'h', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h5.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
   title('P3 Change in Crevasse Height to Width Ratio per Epoch');
   xlabel('Year');
   ylabel('Crevasse H/W');
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;
   % Basal Epoch Change with respect to mid point 
%%
   figure(13)
hold on 
for k = 1:numel(plot_data)
    legend_text =['Crevasse ', num2str(k)];
    size = 70;
    plot(plot_data(k).DX_midpoint, plot_data(k).DH_short_annual,'-o', 'MarkerFaceColor',colors(k,:),...
        'Color',colors(k,:), 'DisplayName', legend_text);
    h = plot(plot_data(k).DX_midpoint, plot_data(k).DH_short_annual, 'o', 'Color','k');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h1 = scatter(plot_data(k).DX_midpoint(1), plot_data(k).DH_short_annual(1), size,'o', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h2 = scatter(plot_data(k).DX_midpoint(2), plot_data(k).DH_short_annual(2), size,'v', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h3 = scatter(plot_data(k).DX_midpoint(3), plot_data(k).DH_short_annual(3), size,'s', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h3.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h4 = scatter(plot_data(k).DX_midpoint(4), plot_data(k).DH_short_annual(4), size,'p', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h4.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
   title('P3 Average Annual Change in Crevasse Height Between Epochs 2010-2011, 2011-2013, 2013-2014, & 2014-2017');
   xlabel('Midpoint distance from Grounding Line (Km)');
   ylabel('Annual \Delta H between Epochs (m/yr)');
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;

   % Change in Crevasse Apex Elevation
   figure(14)
hold on 
for k = 1:numel(plot_data)
    size = 70;
    legend_text =['Crevasse ', num2str(k)];
    plot(plot_data(k).apex_GL_dist, plot_data(k).apex_elev,'-o', 'MarkerFaceColor',colors(k,:),...
        'Color',colors(k,:), 'DisplayName', legend_text);
    h = plot(plot_data(k).apex_GL_dist, plot_data(k).apex_elev, 'o', 'Color','k');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h1 = scatter(plot_data(k).apex_GL_dist(1), plot_data(k).apex_elev(1),size,'o', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h2 = scatter(plot_data(k).apex_GL_dist(2), plot_data(k).apex_elev(2),size, 'v', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h3 = scatter(plot_data(k).apex_GL_dist(3), plot_data(k).apex_elev(3),size, 's', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h3.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h4 = scatter(plot_data(k).apex_GL_dist(4), plot_data(k).apex_elev(4),size, 'p', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h5 = scatter(plot_data(k).apex_GL_dist(5), plot_data(k).apex_elev(5),size, 'h', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h5.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
   title('P3 Crevasse Apex Elevation Change vs. Distance from Grounding Line for 2010, 2011, 2013, 2014, 2017');
   xlabel('Distance from Grounding Line (Km)');
   ylabel('WGS84 Elevation (m)');
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;

   % Average Annual Change in Crevasse Width
figure(15)
hold on 
for k = 1:numel(plot_data)
    size = 70;
    legend_text =['Crevasse ', num2str(k)];
    plot(plot_data(k).DX_midpoint, plot_data(k).DW_annual,'-o', 'MarkerFaceColor',colors(k,:),...
        'Color',colors(k,:), 'DisplayName', legend_text);
    h = plot(plot_data(k).DX_midpoint, plot_data(k).DW_annual, 'o', 'Color','k');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h1 = scatter(plot_data(k).DX_midpoint(1), plot_data(k).DW_annual(1),size,'o', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h2 = scatter(plot_data(k).DX_midpoint(2), plot_data(k).DW_annual(2),size,'v', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h3 = scatter(plot_data(k).DX_midpoint(3), plot_data(k).DW_annual(3),size,'s', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h3.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h4 = scatter(plot_data(k).DX_midpoint(4), plot_data(k).DW_annual(4),size,'p', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h4.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
   title('P3 Average Annual Change in Crevasse Width Between Epochs, 2010-2011, 2011-2013, 2013-2014, & 2014-2017');
   xlabel('Midpoint distance from Grounding Line (Km)');
   ylabel('Annual \Delta W between Epochs (m/yr)');
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;
   
%% Correct
figure(1)
plot(AT_data(1).AT_interp_0km, AT_data(1).Interp_bed);
hold on
plot(AT_data(2).AT_interp_0km,AT_data(2).Interp_bed);
plot(AT_data(3).AT_interp_0km ,AT_data(3).Interp_bed);
plot(AT_data(4).AT_interp_0km ,AT_data(4).Interp_bed);
plot(AT_data(5).AT_interp_0km ,AT_data(5).Interp_bed);
plot(AT_data(1).AT_apex_0km, AT_data(1).crevasse_apex,'o');
plot(AT_data(2).AT_apex_0km, AT_data(2).crevasse_apex,'p');
plot(AT_data(3).AT_apex_0km, AT_data(3).crevasse_apex,'s');
plot(AT_data(4).AT_apex_0km, AT_data(4).crevasse_apex,'x');
plot(AT_data(5).AT_apex_0km, AT_data(5).crevasse_apex,'h');
title('Crevasse Apexes')

figure(2)
plot(AT_data(1).AT_interp_0km, AT_data(1).Interp_bed);
hold on
plot(AT_data(2).AT_interp_0km,AT_data(2).Interp_bed);
plot(AT_data(3).AT_interp_0km ,AT_data(3).Interp_bed);
plot(AT_data(4).AT_interp_0km ,AT_data(4).Interp_bed);
plot(AT_data(5).AT_interp_0km ,AT_data(5).Interp_bed);
plot(AT_data(1).AT_base_0km, AT_data(1).crevasse_base,'o');
plot(AT_data(2).AT_base_0km, AT_data(2).crevasse_base,'p');
plot(AT_data(3).AT_base_0km, AT_data(3).crevasse_base,'s');
plot(AT_data(4).AT_base_0km, AT_data(4).crevasse_base,'x');
plot(AT_data(5).AT_base_0km, AT_data(5).crevasse_base,'h');
title('Crevasse bases')