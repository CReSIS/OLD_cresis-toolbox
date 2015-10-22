% script plot_beam_pattern
%
% Plots the radiation or beam pattern for an antenna configuration
% using the lever_arm.m
%
% Designed for MCoRDS array
% + Should eventually include individual antenna pattern plots rather than
%   assuming isotropic radiators.

param.season_name = '2013_Greenland_P3';
param.radar_name = 'rds';
param.gps_source = 'ATM-final20140301';

physical_constants;
fc = 195e6;
k = 2*pi*fc/c;
tx_weights = hanning(7)';

beam_angles_deg = [-30 0 30];
colors = {'r','g','b'};

%tx_enable_array = {logical([1 1 1 0 0 0 0]), logical([0 0 1 1 1 0 0]), logical([0 0 0 0 1 1 1])};
tx_enable_array = {logical([1 1 1 1 1 1 1])};

normalize = false;

%% Automated Section

figure(1); clf;

h_plot = [];
legend_strs = {};
for wf = 1:length(beam_angles_deg)
  
  if length(tx_enable_array) == 1
    tx_enable = tx_enable_array{1};
  else
    tx_enable = tx_enable_array{wf};
  end
  tx_weights_masked = tx_weights .* tx_enable;
  
  % Get the transmit antenna phase centers
  for tx_chan = 1:7
    tmp_tx_weights = zeros(1,7);
    tmp_tx_weights(tx_chan) = 1;
    phase_centers(:,tx_chan) = lever_arm(param, tmp_tx_weights, tx_chan);
  end
  
  phase_centers = phase_centers - repmat(mean(phase_centers,2), [1 size(phase_centers,2)]);
  
  beam_vec = [0 sind(beam_angles_deg(wf)) cosd(beam_angles_deg(wf))];
  beam_delay = -beam_vec * phase_centers;
  beam_delay = beam_delay / c; % Only include delay relative to nadir
  beam_delay = [beam_delay - mean(beam_delay)];
  
  phase_weights = -beam_delay*fc * 2*pi;
  
  R = 1000;
  phi = linspace(-pi/2, pi/2, 101);
  % phi positive means left looking
  target = [zeros(size(phi)); R*-sin(phi); R*cos(phi)];
  
  beam = zeros(1,size(target,2));
  for idx = 1:size(target,2)
    range = repmat(target(:,idx), [1 size(phase_centers,2)]) - phase_centers;
    range = sqrt(sum(range .* range));
    beam(idx) = sum(tx_weights_masked .* exp(-j*(k*range + phase_weights)));
  end
  
  if normalize
    normalization = max(lp(beam,2));
  else
    normalization = lp(size(phase_centers,2),2);
  end
  
  h_plot(wf) = plot(phi*180/pi, lp(beam,2) - normalization, colors{wf});
  hold on;
  grid on;
  xlabel('Phi (deg)');
  ylabel('Relative power (dB)');
  
  legend_strs{wf} = sprintf('Angle %d', beam_angles_deg(wf));
end
hold off;
legend(h_plot, legend_strs);