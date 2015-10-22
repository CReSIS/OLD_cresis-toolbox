%rds_radiation_pattern_fn = '/mnt/products/tmp/rds_radiation_pattern_2013_Antarctica_Basler_all.mat';
%rds_radiation_pattern_fn = '/mnt/products/tmp/rds_radiation_pattern_2013_Antarctica_Basler_tx.mat';
rds_radiation_pattern_fn = '/cresis/snfs1/scratch/2013_Antarctica_Basler/products/tmp/rds_radiation_pattern_2013_Antarctica_Basler_beams.mat';
tmp = load(rds_radiation_pattern_fn);

rlines = 1:size(tmp.surf_vals{1},2);

color_modes = [0 0 0
    0.3 0.3 0.3
  1 0 0
  1 1 0
  0 1 1
  0 1 0
  0 0 1
  1 0 1];

figure(2); clf;
legend_txt = {};
ant_vals = [];

epri = tmp.epri{1};
for ant_idx = 2:length(tmp.epri)
  epri = intersect(epri,tmp.epri{ant_idx});
end
[epri rlines_master] = intersect(tmp.epri{1}, epri);

fake_shift = zeros(1,length(tmp.surf_vals));
for ant_idx = [1 2 3]%[1 2 3 4 5 6 7 8]%1:length(tmp.surf_vals)
  
  [epri rlines] = intersect(tmp.epri{ant_idx}, epri);
  
  [sort_rolls sort_idxs] = sort(tmp.roll(rlines_master));
  
%   surf_vals_filt = double(abs(tmp.surf_vals{ant_idx}(5,rlines)).^2) / 50;
  surf_vals_filt = double(abs(max(tmp.surf_vals{ant_idx}(:,rlines))).^2) / 50;
  sort_vals = 10*log10(surf_vals_filt) + 30;
  sort_vals = sort_vals(sort_idxs);
  
  ant_vals(:,ant_idx) = medfilt1(sort_vals,301);
%   ant_vals(:,ant_idx) = filter(ones(501,1)/501,1,sort_vals);
  h = plot(sort_rolls*180/pi + fake_shift(ant_idx),ant_vals(:,ant_idx));
  set(h,'Color',color_modes(mod(ant_idx-1,size(color_modes,1))+1,:));
  hold on;
  
  legend_txt{ant_idx} = sprintf('Antenna %d',ant_idx);
  
end
hold off;
xlim([-50 +50])
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


