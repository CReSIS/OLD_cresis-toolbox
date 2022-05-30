%% clip interpolated datasets to grounding line location
for i = 1:numel(XY_data)
    AT_data(i).Interp_bed = XY_data(i).Interp_Bed_Corrected(343:end);
    AT_data(i).Interp_AT_dist = XY_data(i).Interp_AT_Dist(343:end);
    AT_data(i).Interp_AT_dist_0 = AT_data(i).Interp_AT_dist - AT_data(i).Interp_AT_dist(1);
    AT_data(i).InterP_vel_mag = XY_data(i).Interp_Vel_Mag(343:end);
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
AT_data(1).apexes_idx = [4, 12, 14, 17, 20, 27, 31, 38, 50, 62, 86, 120, 165, 222, 273];
%AT_data(1).apex_direction = [r, l, r, r, r, r, l, r, r, l];
% Get index of left and right short sides
AT_data(1).apex_direction_bin_r = [1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0];
AT_data(1).apex_direction_bin_l = [0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1];
AT_data(1).apex_data = AT_data(1).crevasse_apex(AT_data(1).apexes_idx);
AT_data(1).base_data = AT_data(1).crevasse_base(AT_data(1).apexes_idx);
AT_data(1).apex_AT = AT_data(1).AT_dist_apex(AT_data(1).apexes_idx);
AT_data(1).left = [2, 11, 13, 17, 20, 27, 31, 42, 54, 63, 78, 114, 163, 210, 271]; 
AT_data(1).left_base = AT_data(1).AT_dist_base(AT_data(1).left);
AT_data(1).right = [4, 12, 14, 18, 22, 28, 32, 43, 55, 69, 86, 120, 171, 243, 274];
AT_data(1).right_base = AT_data(1).AT_dist_base(AT_data(1).right);


AT_data(2).apexes_idx = [2, 8, 9, 12, 15, 20, 21, 24, 34, 37, 50, 68, 95, 131, 159];
%AT_data(2).apex_direction = [r, r, l, r, l, l, r, r, l, l];
AT_data(2).apex_direction_bin_r = [1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0];
AT_data(2).apex_direction_bin_l = [0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1];
AT_data(2).apex_data = AT_data(2).crevasse_apex(AT_data(2).apexes_idx);
AT_data(2).apex_AT = AT_data(2).AT_dist_apex(AT_data(2).apexes_idx);
AT_data(2).left = [1, 6, 7, 10, 12, 18, 19, 22, 31, 35, 47, 64, 93, 125, 158]; 
AT_data(2).left_base = AT_data(2).AT_dist_base(AT_data(2).left);
AT_data(2).right = [2, 7, 8, 11, 16, 19, 20, 24, 33, 37, 50, 68, 99, 146, 159];
AT_data(2).right_base = AT_data(2).AT_dist_base(AT_data(2).right);
% AT_data(2).base_width = AT_data(2).right_base - AT_data(2).left_base;
% AT_data(2).short_idx_r = find(AT_data(2).apex_direction_bin_r == 1);
% AT_data(2).short_idx_l = find(AT_data(2).apex_direction_bin_l == 1);
% AT_data(2).short_base = [AT_data(2).base_data(AT_data(2).long_idx_r) AT_data(2).base_data(AT_data(2).short_idx_l)];
% AT_data(2).long_idx_r = find(AT_data(2).apex_direction_bin_r == 0);
% AT_data(2).long_idx_l = find(AT_data(2).apex_direction_bin_l == 0);
% AT_data(2).long_base = [AT_data(2).base_data(AT_data(2).long_idx_r) AT_data(2).base_data(AT_data(2).long_idx_l)];
% AT_data(1).DH_long = AT_data(1).apex_data - AT_data(1).long_base;
% AT_data(1).DH_short = AT_data(1).apex_data - AT_data(1).short_base;

AT_data(3).apexes_idx = [1, 3, 3, 11, 11, 19, 21, 28, 35, 47, 61, 82, 113, 175 190];
%AT_data(3).apex_direction = [r, r, r, r, r, r, r, l, r, r];
AT_data(3).apex_direction_bin_r = [1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0];
AT_data(3).apex_direction_bin_l = [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1];
AT_data(3).apex_data = AT_data(3).crevasse_apex(AT_data(3).apexes_idx);
AT_data(3).apex_AT = AT_data(3).AT_dist_apex(AT_data(3).apexes_idx);
AT_data(3).left_base_1 = [AT_data(3).AT_dist_base(1)];  
AT_data(3).left = [2, 2, 10, 10, 18, 21, 25, 35, 42, 59, 74, 110, 147, 186]; 
AT_data(3).left_base = [AT_data(3).left_base_1 AT_data(3).AT_dist_base(AT_data(3).left)];
AT_data(3).left_base(13) = AT_data(3).left_base(13) + (0.2876*(1e3));
AT_data(3).right = [2, 6, 6, 16, 16, 21, 23, 31, 42, 51, 63, 82, 132, 178, 190];
AT_data(3).right_base = AT_data(3).AT_dist_base(AT_data(3).right);
% AT_data(3).base_width = AT_data(3).right_base - AT_data(3).left_base;
% AT_data(3).short_idx_r = find(AT_data(3).apex_direction_bin_r == 1);
% AT_data(3).short_idx_l = find(AT_data(3).apex_direction_bin_l == 1);
% AT_data(3).short_base = [AT_data(3).base_data(AT_data(3).long_idx_r) AT_data(3).base_data(AT_data(3).short_idx_l)];
% AT_data(3).long_idx_r = find(AT_data(3).apex_direction_bin_r == 0);
% AT_data(3).long_idx_l = find(AT_data(3).apex_direction_bin_l == 0);
% AT_data(3).long_base = [AT_data(3).base_data(AT_data(3).long_idx_r) AT_data(3).base_data(AT_data(3).long_idx_l)];
% AT_data(3).DH_long = AT_data(3).apex_data - AT_data(3).long_base;
% AT_data(3).DH_short = AT_data(3).apex_data - AT_data(3).short_base;

for i = 1:3
    AT_data(i).apex_data = AT_data(i).crevasse_apex(AT_data(i).apexes_idx);
    AT_data(i).base_data = AT_data(i).crevasse_base(AT_data(i).apexes_idx);
    AT_data(i).base_width = AT_data(i).right_base - AT_data(i).left_base;
    AT_data(i).short_idx_r = find(AT_data(i).apex_direction_bin_r == 1);
    AT_data(i).short_idx_l = find(AT_data(i).apex_direction_bin_l == 1);
    AT_data(i).short_base = [AT_data(i).base_data(AT_data(i).short_idx_r) AT_data(i).base_data(AT_data(i).short_idx_l)];
    AT_data(i).long_idx_r = find(AT_data(i).apex_direction_bin_r == 0);
    AT_data(i).long_idx_l = find(AT_data(i).apex_direction_bin_l == 0);
    AT_data(i).long_base = [AT_data(i).base_data(AT_data(i).long_idx_r) AT_data(i).base_data(AT_data(i).long_idx_l)];
    AT_data(i).DH_long = AT_data(i).apex_data - AT_data(i).long_base;
    AT_data(i).DH_short = AT_data(i).apex_data - AT_data(i).short_base;
end

% Define time period off summations for annual correction
time_offset.Period_1 = time_offset.P13 + time_offset.P14_1;
time_offset.Period_2 = time_offset.P14_2 + time_offset.P15 + time_offset.P16 + time_offset.P17;

% loop for taking manual data and splitting into groupings per crevasse
for i = 1:length(AT_data(3).apexes_idx)
    plot_data(i).apex_elev = [AT_data(1).apex_data(i), AT_data(3).apex_data(i), AT_data(3).apex_data(i)];
    plot_data(i).apex_GL_dist = [AT_data(1).apex_AT(i)/1e3, AT_data(2).apex_AT(i)/1e3, AT_data(3).apex_AT(i)/1e3];
    plot_data(i).DH_short = [AT_data(1).DH_short(i), AT_data(2).DH_short(i), AT_data(3).DH_short(i)];
    plot_data(i).DH_long = [AT_data(1).DH_long(i), AT_data(2).DH_long(i), AT_data(3).DH_long(i)];
    plot_data(i).widths = [AT_data(1).base_width(i), AT_data(2).base_width(i), AT_data(3).base_width(i)];
    plot_data(i).DW_annual = [((AT_data(2).base_width(i) - AT_data(1).base_width(i))/(time_offset.Period_1)), ...
        ((AT_data(3).base_width(i) - AT_data(2).base_width(i))/(time_offset.Period_2))];
    plot_data(i).DH_short_annual = [((AT_data(2).DH_short(i) - AT_data(1).DH_short(i))/(time_offset.Period_1)), ...
        ((AT_data(3).DH_short(i) - AT_data(2).DH_short(i))/(time_offset.Period_2))];
    plot_data(i).DH_long_annual = [((AT_data(2).DH_long(i) - AT_data(1).DH_long(i))/(time_offset.Period_1)), ...
        ((AT_data(3).DH_long(i) - AT_data(2).DH_long(i))/(time_offset.Period_2))];
    plot_data(i).DX = [(AT_data(2).apex_AT(i)/1e3 - AT_data(1).apex_AT(i)/1e3), (AT_data(3).apex_AT(i)/1e3 - AT_data(2).apex_AT(i)/1e3)];
    plot_data(i).DX_yr = [((AT_data(2).apex_AT(i)/1e3 - AT_data(1).apex_AT(i)/1e3)/time_offset.Period_1),...
        ((AT_data(3).apex_AT(i)/1e3 - AT_data(2).apex_AT(i)/1e3)/time_offset.Period_2)];
    plot_data(i).DX_midpoint = [((AT_data(2).apex_AT(i)/1e3 + AT_data(1).apex_AT(i)/1e3)/2), ...
        ((AT_data(3).apex_AT(i)/1e3 + AT_data(2).apex_AT(i)/1e3)/2)];
end
%% Figures for Poster
   % Basal Epoch Change with respect to mid point 
colors = colormap(flipud(turbo(length(plot_data))));
figure(10)
hold on 
for k = 1:numel(plot_data)
    legend_text =['Crevasse ', num2str(k)];
    sizes1 = plot_data(i).DW_annual(1);
    sizes2 = plot_data(i).DW_annual(2);
    plot(plot_data(k).DX_midpoint, plot_data(k).DH_short_annual,'-o', 'MarkerFaceColor',colors(k,:),...
        'Color',colors(k,:), 'DisplayName', legend_text);
    h = plot(plot_data(k).DX_midpoint, plot_data(k).DH_short_annual, 'o', 'Color','k');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h1 = scatter(plot_data(k).DX_midpoint(1), plot_data(k).DH_short_annual(1), sizes1,'o', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h2 = scatter(plot_data(k).DX_midpoint(2), plot_data(k).DH_short_annual(2), sizes2,'v', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
   title('P2 Average Annual Change in Crevasse Height Between Epochs 2013-2014 & 2014-2017');
   xlabel('Midpoint distance from Grounding Line (Km)');
   ylabel('Annual \Delta H between Epochs (m/yr)');
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;

   % Change in Crevasse Apex Elevation
   figure(11)
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
end
   title('P2 Crevasse Apex Elevation Change vs. Distance from Grounding Line for 2013, 2014, 2017');
   xlabel('Distance from Grounding Line (Km)');
   ylabel('WGS84 Elevation (m)');
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;

   % Average Annual Change in Crevasse Width
figure(12)
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
end
   title('P2 Average Annual Change in Crevasse Width Between Epochs, 2013-2014 & 2014-2017');
   xlabel('Midpoint distance from Grounding Line (Km)');
   ylabel('Annual \Delta W between Epochs (m/yr)');
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;
%%
% define color ramp to use
% total crevasse change per year with respect to grounding line
colors = colormap(flipud(turbo(length(plot_data))));
figure(21)
hold on 
for k = 1:numel(plot_data)
    legend_text =['Crevasse ', num2str(k)];
    sizes = plot_data(i).widths;
    plot(plot_data(k).apex_GL_dist, plot_data(k).DH_short,'-o', 'MarkerFaceColor',colors(k,:),...
        'Color',colors(k,:), 'DisplayName', legend_text);
    h = plot(plot_data(k).apex_GL_dist, plot_data(k).DH_short, 'o', 'Color','k');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h2 = scatter(plot_data(k).apex_GL_dist, plot_data(k).DH_short, sizes, 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
   title('P2 Crevasse Height and Width Change as function of Distance from Grounding Line ');
   xlabel('Distance from Grounding Line (Km)');
   ylabel('\Delta H over observered period (m)');
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;
% Basal Epoch Change with respect to velocity   
figure(22)
hold on 
for k = 1:numel(plot_data)
    legend_text =['Crevasse ', num2str(k)];
    sizes1 = plot_data(i).DW_annual(1);
    sizes2 = plot_data(i).DW_annual(2);
    plot(plot_data(k).DX_yr, plot_data(k).DH_short_annual,'-o', 'MarkerFaceColor',colors(k,:),...
        'Color',colors(k,:), 'DisplayName', legend_text);
    h = plot(plot_data(k).DX_yr, plot_data(k).DH_short_annual, 'o', 'Color','k');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h1 = scatter(plot_data(k).DX_yr(1), plot_data(k).DH_short_annual(1), sizes1,'o', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h2 = scatter(plot_data(k).DX_yr(2), plot_data(k).DH_short_annual(2), sizes2,'v', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
   title('P2 Change in Crevasse Height and Width Between Epochs');
   xlabel('Average Annual Velocity (Km/yr)');
   ylabel('Annual \Delta H between Epochs (m/yr)');
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;  

   % Basal Epoch Change with respect to mid point 
   figure(29)
hold on 
for k = 1:numel(plot_data)
    legend_text =['Crevasse ', num2str(k)];
    sizes1 = plot_data(i).DW_annual(1);
    sizes2 = plot_data(i).DW_annual(1);
    plot(plot_data(k).DX_midpoint, plot_data(k).DH_short_annual,'-o', 'MarkerFaceColor',colors(k,:),...
        'Color',colors(k,:), 'DisplayName', legend_text);
    h = plot(plot_data(k).DX_midpoint, plot_data(k).DH_short_annual, 'o', 'Color','k');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h1 = scatter(plot_data(k).DX_midpoint(1), plot_data(k).DH_short_annual(1), sizes1,'o', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h2 = scatter(plot_data(k).DX_midpoint(2), plot_data(k).DH_short_annual(2), sizes2,'v', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
   title('P2 Change in Crevasse Height and Width Between Epochs');
   xlabel('Average Annual Velocity (Km/yr)');
   ylabel('Annual \Delta H between Epochs (m/yr)');
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show; 
   
   % Normalized Basal Crevasse change (not annually normalized though)
   total_std_heights = std([plot_data.DH_short]);
   total_std_widths = std([plot_data.widths]);
   mean_height = mean([plot_data.DH_short]);
   mean_width = mean([plot_data.widths]);
   figure(23)
   hold on 
for k = 1:numel(plot_data)
    legend_text =['Crevasse ', num2str(k)];
    sizes = abs(((plot_data(i).widths - mean_width)/total_std_widths)*100);
    plot(plot_data(k).apex_GL_dist, ((plot_data(k).DH_short - mean_height)/total_std_heights),'-o', 'MarkerFaceColor',colors(k,:),...
        'Color',colors(k,:), 'DisplayName', legend_text);
    h = plot(plot_data(k).apex_GL_dist, ((plot_data(k).DH_short - mean_height)/total_std_heights), 'o', 'Color','k');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h2 = scatter(plot_data(k).apex_GL_dist, ((plot_data(k).DH_short - mean_height)/total_std_heights), sizes, 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
   title('P2 Normalized Crevasse Height and Width Change as function of Distance from Grounding Line ');
   xlabel('Distance from Grounding Line (Km)');
   ylabel('\Delta H over observered period (m)');
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;
   
   % Annual Width vs. Height
   figure(24)
   hold on 
for k = 1:numel(plot_data)
    legend_text =['Crevasse ', num2str(k)];
    sizes = plot_data(i).widths;
    plot(plot_data(k).DW_annual, plot_data(k).DH_short_annual,'-o', 'MarkerFaceColor',colors(k,:),...
        'Color',colors(k,:), 'DisplayName', legend_text);
    h = plot(plot_data(k).DW_annual, plot_data(k).DH_short_annual, 'o', 'Color','k');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h1 = scatter(plot_data(k).DW_annual(1), plot_data(k).DH_short_annual(1), sizes1,'o', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h2 = scatter(plot_data(k).DW_annual(2), plot_data(k).DH_short_annual(2), sizes2,'v', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
   title('P2 Crevasse Height vs. Width Change ');
   xlabel('Annual \Delta W (m/yr)');
   ylabel('Annual \Delta H (m/yr)');
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;
   
%%  Basal Crevasse Apex Postion change 
   figure(25)
hold on 
for k = 1:numel(plot_data)
    legend_text =['Crevasse ', num2str(k)];
    plot(plot_data(k).apex_GL_dist, plot_data(k).apex_elev,'-o', 'MarkerFaceColor',colors(k,:),...
        'Color',colors(k,:), 'DisplayName', legend_text);
    h = plot(plot_data(k).apex_GL_dist, plot_data(k).apex_elev, 'o', 'Color','k');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h1 = scatter(plot_data(k).apex_GL_dist(1), plot_data(k).apex_elev(1), 'o', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h2 = scatter(plot_data(k).apex_GL_dist(2), plot_data(k).apex_elev(2), 'o', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h3 = scatter(plot_data(k).apex_GL_dist(3), plot_data(k).apex_elev(3), 'o', 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    h3.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
   title('P2 Crevasse Apex Position as function of Distance from Grounding Line ');
   xlabel('Distance from Grounding Line (Km)');
   ylabel('WGS84 Elevation (m)');
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;
%% 
   colors = colormap(flipud(turbo(length(DH_Max_points))));
   total_std = std([DH_width_points.points]);
   total_mean = mean([DH_width_points.points]);
   figure(20)
   hold on
   for k = 1:numel(DH_pass)
       text_1 = ['Crevasse ', num2str(k)];
       %sizes = (DH_width_points(k).points + 1)*100;
       sizes = abs(((DH_width_points(k).points - total_mean)/total_std)*100);
       plot(DH_pass(k).points, DH_points(k).points, '-o', 'MarkerFaceColor',colors(k,:), 'Color',colors(k,:),  'DisplayName', text_1);
       h = plot(DH_pass(k).points, DH_points(k).points, 'o', 'Color','k');
       h.Annotation.LegendInformation.IconDisplayStyle = 'off';
       h2 = scatter(DH_pass(k).points, DH_points(k).points, sizes, 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
       h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
   end
   title('Crevasse Apex Elevation change with Respect to WGS84 0-surface as a function of Distance from Grounding Line and Crevasse Width');
   xlabel('Distance from Grounding Line (Km)');
   ylabel('Apex Elevation (m)');
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;
 
   %%
figure(1)
plot(AT_data(1).Interp_AT_vel - AT_data(1).Interp_AT_vel(1) ,AT_data(1).Interp_bed);
hold on
plot(AT_data(2).Interp_AT_vel - AT_data(2).Interp_AT_vel(1) ,AT_data(2).Interp_bed);
plot(AT_data(3).Interp_AT_vel - AT_data(3).Interp_AT_vel(1) ,AT_data(3).Interp_bed);
plot(AT_data(1).AT_dist_apex_vel - AT_data(1).Interp_AT_vel(1), AT_data(1).crevasse_apex,'o');
plot(AT_data(2).AT_dist_apex_vel - AT_data(2).Interp_AT_vel(1), AT_data(2).crevasse_apex,'o');
plot(AT_data(3).AT_dist_apex_vel - AT_data(3).Interp_AT_vel(1), AT_data(3).crevasse_apex,'o');
   
%% Correct
figure(1999)
plot(AT_data(1).Interp_AT_vel/1e3 - AT_data(1).Interp_AT_vel(1)/1e3, AT_data(1).Interp_bed)
hold on
plot(AT_data(2).Interp_AT_vel/1e3 - AT_data(2).Interp_AT_vel(1)/1e3, AT_data(2).Interp_bed)
plot(AT_data(3).Interp_AT_vel/1e3 - AT_data(3).Interp_AT_vel(1)/1e3, AT_data(3).Interp_bed)
plot(AT_data(1).AT_dist_apex/1e3 - AT_data(1).Interp_AT_vel(1)/1e3, AT_data(1).crevasse_apex,'o')
plot(AT_data(2).AT_dist_apex/1e3 - AT_data(2).Interp_AT_vel(1)/1e3, AT_data(2).crevasse_apex,'d')
plot(AT_data(3).AT_dist_apex/1e3 - AT_data(3).Interp_AT_vel(1)/1e3, AT_data(3).crevasse_apex,'x')
title('Crevasse Apexes')

figure(2000)
plot(AT_data(1).Interp_AT_vel/1e3 - AT_data(1).Interp_AT_vel(1)/1e3, AT_data(1).Interp_bed)
hold on
plot(AT_data(2).Interp_AT_vel/1e3 - AT_data(2).Interp_AT_vel(2)/1e3, AT_data(2).Interp_bed)
plot(AT_data(3).Interp_AT_vel/1e3 - AT_data(3).Interp_AT_vel(3)/1e3, AT_data(3).Interp_bed)
plot(AT_data(1).AT_dist_base/1e3 - AT_data(1).Interp_AT_vel(1)/1e3, AT_data(1).crevasse_base,'o')
plot(AT_data(2).AT_dist_base/1e3 - AT_data(2).Interp_AT_vel(1)/1e3, AT_data(2).crevasse_base,'d')
plot(AT_data(3).AT_dist_base/1e3 - AT_data(3).Interp_AT_vel(1)/1e3, AT_data(3).crevasse_base,'x')
title('Crevasse Apexes')

figure(2021)
plot(AT_data(1).Interp_AT_vel/1e3 - AT_data(1).Interp_AT_vel(1)/1e3, AT_data(1).Interp_bed)
hold on
plot(AT_data(2).Interp_AT_vel/1e3 - AT_data(2).Interp_AT_vel(1)/1e3, AT_data(2).Interp_bed)
plot(AT_data(3).Interp_AT_vel/1e3 - AT_data(3).Interp_AT_vel(1)/1e3, AT_data(3).Interp_bed)
plot(AT_data(1).AT_dist_apex/1e3 - AT_data(1).Interp_AT_vel(1)/1e3, AT_data(1).crevasse_apex,'o')
plot(AT_data(2).AT_dist_apex/1e3 - AT_data(2).Interp_AT_vel(1)/1e3, AT_data(2).crevasse_apex,'d')
plot(AT_data(3).AT_dist_apex/1e3 - AT_data(3).Interp_AT_vel(1)/1e3, AT_data(3).crevasse_apex,'x')
title('Crevasse Apexes')
%%
for i =1:3
    AT_data(i).Aligned_at_0 = AT_data(i).AT_dist_base_vel/1e3 - AT_data(i).Interp_AT_vel(1)/1e3;
end
%%
figure(2002)
plot(AT_data(1).Interp_AT_dist, AT_data(1).Interp_bed)
hold on
plot(AT_data(2).Interp_AT_dist, AT_data(2).Interp_bed)
plot(AT_data(3).Interp_AT_dist, AT_data(3).Interp_bed)
% plot(AT_data(1).AT_dist_apex, AT_data(1).crevasse_apex,'o')
% plot(AT_data(2).AT_dist_apex, AT_data(2).crevasse_apex,'d')
% plot(AT_data(1).Interp_AT_dist, AT_data(3).Interp_bed)
ylim([-500 0]);

figure(2001)
plot(AT_data(1).Interp_AT_dist, AT_data(1).Interp_bed)
hold on
plot(AT_data(2).Interp_AT_dist, AT_data(2).Interp_bed)
plot(AT_data(3).Interp_AT_dist, AT_data(3).Interp_bed)
ylim([-500 0]);
%% AFTER MANUAL EDITTING EXPORT AND SAVE AS CSV
% after mannually editting of each profile, save as new variable and export
% variable as a csv 
% DO NO DO IF CSV ALREADY EXISTS
% COMMENT SECTION OUT AFTER FINALS PICKS AS TO NO OVERWRITE EXISTING FILES
crevasse_2013 = AT_data.crevasse.base_pic_P13;
crevasse_2014 = AT_data.crevasse.base_pic_P14;
crevasse_2017 = AT_data.crevasse.base_pic_P17;
writematrix(crevasse_2013, 'P1E_Crevasse_2013.csv','Delimiter','comma');
writematrix(crevasse_2014, 'P1E_Crevasse_2014.csv','Delimiter','comma');
writematrix(crevasse_2017, 'P1E_Crevasse_2017.csv','Delimiter','comma');
%% Reload new csv in after editting
AT_data.crevasse.base_pic_P11 = csvread('P1E_Crevasse_2011.csv');
AT_data.crevasse.base_pic_P14 = csvread('P1E_Crevasse_2014.csv');
AT_data.crevasse.base_pic_P15 = csvread('P1E_Crevasse_2015.csv');
AT_data.crevasse.base_pic_P17 = csvread('P1E_Crevasse_2017.csv');
AT_data.crevasse.base_pic_P18 = csvread('P1E_Crevasse_2018.csv');
AT_data.crevasse.base_pic_P19 = csvread('P1E_Crevasse_2019.csv');
%% 2011 Basal Crevasse picks
P1_index = []
%% 2014 Basal Crevasse picks

%% 2015 Basal Crevasse picks

%% 2017 Basal Crevasse picks

%% 2018 Basal Crevasse Picks

%% 2019 Basal Crevasse picks







%% Basal Crevasse Plotting
figure(10)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P11,'r');
hold on;
h2 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P14,'b');
h3 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P18,'g');
h4 = plot(AT_data.query_array.P11/1e3, AT_data.crevasse.apex_pic_P11,'k*');
h5 = plot(AT_data.query_array.P11/1e3, AT_data.crevasse.apex_pic_P14,'k*');
h6 = plot(AT_data.query_array.P11/1e3, AT_data.crevasse.apex_pic_P18,'k*');
title('Basal Crevasse Apex Picks - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('corrected 2011', 'corrected 2014', 'corrected 2018', ...
  'Apex Pick 2011', 'Apex Pick 2014', 'Apex Pick 2018', ...
  'Location', 'southeast');
hold off

figure(11)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P11,'r');
hold on;
h2 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P14,'b');
h3 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P18,'g');
h4 = plot(AT_data.query_array.P11/1e3, AT_data.crevasse.base_pic_P11,'k*');
h5 = plot(AT_data.query_array.P11/1e3, AT_data.crevasse.base_pic_P14,'k*');
h6 = plot(AT_data.query_array.P11/1e3, AT_data.crevasse.base_pic_P18,'k*');
title('Basal Crevasse Base Picks - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('corrected 2011', 'corrected 2014', 'corrected 2018', ...
  'Base Pick 2011', 'Base Pick 2014', 'Base Pick 2018', ...
  'Location', 'southeast');
hold off;

figure(12)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P11,'r');
hold on
h2 = plot(AT_data.query_array.P11/1e3, AT_data.crevasse.base_pic_P11,'k*');
title('Basal Crevasse Base Picks (2011) - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('corrected 2011', 'Base Pick 2011', 'Location', 'southeast');

figure(13)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P14,'b');
hold on
h2 = plot(AT_data.query_array.P11/1e3, AT_data.crevasse.base_pic_P14,'k*');
title('Basal Crevasse Base Picks (2014) - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('corrected 2014', 'Base Pick 2014', 'Location', 'southeast');

figure(14)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P18,'g');
hold on
h2 = plot(AT_data.query_array.P11/1e3, AT_data.crevasse.base_pic_P18,'k*');
title('Basal Crevasse Base Picks (2014) - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('corrected 2018', 'Base Pick 2018', 'Location', 'southeast');

crevasse1 = AT_data.interp_data.P11(479:497);
along_track_crevasse1 = AT_data.query_array.P11(479:497)/1e3;
cross_area = cumtrapz(along_track_crevasse1, crevasse1);
figure(15)
h1 = plot(along_track_crevasse1, crevasse1);
cross_area1 = patch(along_track_crevasse1, crevasse1,'b','LineWidth',1.5);
title('Basal Crevasse Base Pick 1 (2011) - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('basal crevasse 1', 'Location', 'southeast');

