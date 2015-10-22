
if 1
  tmp = load('/cresis/scratch1/paden/roll_measurements_wf4.mat');
  tmp2 = load('/cresis/scratch1/paden/roll_measurements.mat');
  
  % Fix 
  tmp2.surf_vals{4} = tmp.surf_vals{1};
else
  tmp2 = load('/cresis/scratch1/paden/roll_measurements_ant1.mat');
end

% Switch 2 and 5
switch_tmp = tmp2.surf_vals{2};
tmp2.surf_vals{2} = tmp2.surf_vals{5};
tmp2.surf_vals{5} = switch_tmp;

if 0
  figure(1); clf;
  plot(lp(medfilt1(double(tmp2.surf_vals{1}),51)),'g')
  hold on;
  plot(lp(medfilt1(double(tmp2.surf_vals{2}),51)),'r')
  plot(lp(medfilt1(double(tmp2.surf_vals{7}),51)),'k')
  plot(lp(medfilt1(double(tmp2.surf_vals{4}),51)),'b')
  hold off;
end

rlines = 1:133500;

if 0
  surf_vals_filt = medfilt1(double(tmp2.surf_vals{ant_idx}(rlines)) / 50,101);
  
  figure(1); clf;
  plot(roll(rlines)*180/pi, 10*log10(surf_vals_filt) + 30,'.');
  grid on;
  xlabel('Roll (deg)');
  ylabel('Pr (dBm)');
end

color_modes = [0 0 0
%   0.3 0.3 0.3
  1 0 0
  1 1 0
  0 1 1
  0 1 0
  0 0 1
  1 0 1];


figure(1); clf;
legend_txt = {};
ant_vals = [];
for ant_idx = 1:length(tmp2.surf_vals)
  
  [sort_rolls sort_idxs] = sort(tmp2.roll(rlines));
  
  surf_vals_filt = double(tmp2.surf_vals{ant_idx}(rlines)) / 50;
  sort_vals = 10*log10(surf_vals_filt(rlines)) + 30;
  sort_vals = sort_vals(sort_idxs);
  
  ant_vals(:,ant_idx) = medfilt1(sort_vals,301);
  h = plot(sort_rolls*180/pi,ant_vals(:,ant_idx));
  set(h,'Color',color_modes(ant_idx,:));
  hold on;
  
  legend_txt{ant_idx} = sprintf('Antenna %d',ant_idx);
  
  
end
hold off;
xlim([-37 +37])
legend(legend_txt);
grid on;

roll_angle = sort_rolls;
save('/cresis/scratch1/paden/radiation_pattern.mat','ant_vals','roll_angle');




