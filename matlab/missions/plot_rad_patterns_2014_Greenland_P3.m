% script plot_rad_patterns_2014_Greenland_P3
%
% Uses output from basic_rds_radiation_pattern_*.m functions
%
% Author: John Paden

output_fn = '/cresis/scratch1/paden/roll_measurements_2014_Greenland_P3.mat';
tmp = load(output_fn);

rlines = 1:size(tmp.surf_vals{1},2);

color_modes = [0 0 0
    %0.3 0.3 0.3
  1 0 0
  1 1 0
  0 1 1
  0 1 0
  0 0 1
  1 0 1];

figure(2); clf;
legend_txt = {};
ant_vals = [];

% epri = tmp.epri{1};
% for ant_idx = 2:length(tmp.epri)
%   epri = intersect(epri,tmp.epri{ant_idx});
% end
% [epri rlines_master] = intersect(tmp.epri{1}, epri);
rlines_master = rlines;

fake_shift = zeros(1,length(tmp.surf_vals));

antennas = [2 3 1 4 5 6 7];
antennas = [8];
for ant_idx = 1:length(antennas)%1:length(tmp.surf_vals)
  
  meas_idx = antennas(ant_idx);
%   [epri rlines] = intersect(tmp.epri{ant_idx}, epri);
  
  [sort_rolls sort_idxs] = sort(tmp.roll(rlines_master));
  
  surf_vals_filt = double(abs(tmp.surf_vals{meas_idx}(5,rlines)).^2) / 50;
  sort_vals = 10*log10(surf_vals_filt) + 30;
  sort_vals = sort_vals(sort_idxs);
  
  ant_vals(:,meas_idx) = medfilt1(sort_vals,301);
  h = plot(sort_rolls*180/pi + fake_shift(meas_idx), ant_vals(:,meas_idx) );
  set(h,'Color',color_modes(mod(ant_idx-1,size(color_modes,1))+1,:));
  hold on;
  
  legend_txt{ant_idx} = sprintf('Antenna %d',ant_idx);
  
end
hold off;
xlim([-35 +35])
legend(legend_txt);
grid on;
xlabel('Roll (deg)')
ylabel('Relative power (dB)')

roll_angle = sort_rolls;
% save('/cresis/scratch1/paden/radiation_pattern.mat','ant_vals','roll_angle');

return;











% Form interferogram (couple options)
complex_data = filter2(ones(5)/25,tmp.surf_vals{2} .* conj(tmp.surf_vals{1}));
% complex_data = data{1}(:,:,1) .* conj(data{1}(:,:,2))
% Plot interferogram
imagesc(hsv_plot(complex_data,-70));
colormap(hsv(256))
h_colorbar = colorbar;
caxis([-pi pi])
set(get(h_colorbar,'ylabel'),'string','angle (rad)')


