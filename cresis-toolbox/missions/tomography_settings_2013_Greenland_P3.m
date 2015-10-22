
% 7-element array
beamwidth_in_deg = 1/3.5 * 180/pi

steering_angle = beamwidth_in_deg/2

steering_angle = [15 -22.5];

phase_settings_from_GUI = [64	38	0	67	247	32	102];

param = [];
param.season_name = '2013_Greenland_P3';
param.radar_name = 'mcords3';
param.gps_source = 'ATM-final';
tx_weights = [1 1 1 1 1 1 1];
rxchannels = [1 2 3 4 5 6 7];

for rx = 1:length(rxchannels)
  [phase_center(:,rx)] = lever_arm(param, tx_weights, rxchannels(rx));
end

fc = 195e6;
% [theta,sv] = array_proc_sv({'theta',[-15 0 17.5]/180*pi}, fc, phase_center(2,:).', phase_center(3,:).');
[theta,sv] = array_proc_sv({'theta',[-15 0 10]/180*pi}, fc, phase_center(2,:).', phase_center(3,:).');

angle_rad = angle(sv);
angle_rad(:,2:3) = unwrap(angle_rad(:,2:3));
angle_rad(:,1) = flipud(unwrap(flipud(angle_rad(:,1))));
angle_rad = angle_rad - angle_rad(4,2);
angle_deg = angle_rad*180/pi

angle_deg_left = -(angle_deg(:,1) - angle_deg(:,2)).'
new_GUI_phase_settings = round(phase_settings_from_GUI + angle_deg_left)
angle_deg_right = -(angle_deg(:,3) - angle_deg(:,2)).'
new_GUI_phase_settings = round(phase_settings_from_GUI + angle_deg_right)

