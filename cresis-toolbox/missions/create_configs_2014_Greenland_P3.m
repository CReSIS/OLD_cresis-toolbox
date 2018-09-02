% script create_configs_2014_Greenland_P3
%
% Creates NI radar depth sounder settings

% Define waveforms
base_dir = '/cresis/snfs1/scratch/paden/waveforms/';
% base_dir = 'c:\tmp\';
final_DDS_phase = [179.7	-174	147.4	0	110.5	133.7	116.5 0];
final_DDS_amp = [11849	14598	17772	30000	16162	13620	13140	0];
final_DDS_time = [3.7	3.35	3.67	-0.63	-0.45	-1.38	3 0];
f0 = 180e6;
f1 = 210e6;

Hwindow_orig = [1 1 1 1 1 1 1 1]; % Desired window created during transmit calibration
Hwindow35 = [0.2008    0.5196    0.8557    1.0000    0.8557    0.5196    0.2008 0]; % chebwin(7,35)

physical_constants;
param.season_name = '2014_Greenland_P3';
param.radar_name = 'rds';
param.gps_source = 'atm-final20140301';
clear phase_centers;
for tx_chan = 1:7
  tx_weights = zeros(1,7);
  tx_weights(tx_chan) = 1;
  phase_centers(:,tx_chan) = lever_arm(param, tx_weights, tx_chan);
end

%% Survey Mode
% <1500 m, 500-2000 ft AGL
param = struct('radar_name','mcords3','num_chan',16,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0'); 
param.max_tx = 30000; param.max_data_rate = 135; param.flight_hours = 7; param.sys_delay = 12.24e-6;
param.max_duty_cycle = 0.2;
param.create_IQ = false;
param.tg.altitude_guard = 1000*12*2.54/100;
param.tg.staged_recording = true;
param.tg.Haltitude = 1000*12*2.54/100;
param.tg.Hice_thick = 2500;
param.fn = fullfile(base_dir,'survey_mode_1us_2wf_1500mthick.xml');
param.prf = 12000;
param.presums = [3 33];
param.wfs(1).atten = 30;
param.wfs(2).atten = 0;
DDS_amp = final_DDS_amp;
DDS_amp = round(DDS_amp * param.max_tx / max(DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.2;
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 1e-6;
param.phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
[param.wfs(:).tx_mask] = deal([0 0 0 0 0 0 0 0]);
write_cresis_xml(param);

%% Survey Mode
% <1500 m, 500-2000 ft AGL
param = struct('radar_name','mcords3','num_chan',16,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0'); 
param.max_tx = 30000; param.max_data_rate = 135; param.flight_hours = 7; param.sys_delay = 12.24e-6;
param.max_duty_cycle = 0.2;
param.create_IQ = false;
param.tg.altitude_guard = 1000*12*2.54/100;
param.tg.staged_recording = true;
param.tg.Haltitude = 1000*12*2.54/100;
param.tg.Hice_thick = 2500;
param.fn = fullfile(base_dir,'survey_mode_3us_2wf_1500mthick.xml');
param.prf = 12000;
param.presums = [15 21];
param.wfs(1).atten = 30;
param.wfs(2).atten = 0;
DDS_amp = final_DDS_amp;
DDS_amp = round(DDS_amp * param.max_tx / max(DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.2;
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 3e-6;
param.phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
[param.wfs(:).tx_mask] = deal([0 0 0 0 0 0 0 0]);
write_cresis_xml(param);

%% Survey Mode
% <3500 m, 500-2000 ft AGL
param = struct('radar_name','mcords3','num_chan',16,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0');  
param.max_tx = 30000; param.max_data_rate = 135; param.flight_hours = 7; param.sys_delay = 12.24e-6;
param.max_duty_cycle = 0.2;
param.create_IQ = false;
param.tg.altitude_guard = 1000*12*2.54/100;
param.tg.staged_recording = true;
param.tg.Haltitude = 1000*12*2.54/100;
param.tg.Hice_thick = 3500;
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_3000mthick.xml');
param.prf = 12000;
param.presums = [3 3 29];
param.wfs(1).atten = 30;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
DDS_amp = final_DDS_amp;
DDS_amp = round(DDS_amp * param.max_tx / max(DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.1;
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 3e-6;
param.wfs(3).Tpd = 10e-6;
param.phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
[param.wfs(:).tx_mask] = deal([0 0 0 0 0 0 0 0]);
write_cresis_xml(param);

%% Survey Mode High Altitude
% <3500 m, 10000 ft AGL
param = struct('radar_name','mcords3','num_chan',16,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0');  
param.max_tx = 30000; param.max_data_rate = 135; param.flight_hours = 7; param.sys_delay = 12.24e-6;
param.max_duty_cycle = 0.2;
param.create_IQ = false;
param.tg.altitude_guard = 1500*12*2.54/100;
param.tg.staged_recording = true;
param.tg.Haltitude = 10000*12*2.54/100;
param.tg.Hice_thick = 3500;
param.fn = fullfile(base_dir,'survey_mode_10us_2wf_3000mthick_high_altitude.xml');
param.prf = 10000;
param.presums = [3 29];
param.wfs(1).atten = 10;
param.wfs(2).atten = 0;
DDS_amp = final_DDS_amp;
DDS_amp = round(DDS_amp * param.max_tx / max(DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.1;
param.wfs(1).Tpd = 3e-6;
param.wfs(2).Tpd = 10e-6;
param.phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
[param.wfs(:).tx_mask] = deal([0 0 0 0 0 0 0 0]);
write_cresis_xml(param);

%% Noise Mode
% Receivers should be terminated in 50 Ohms
param = struct('radar_name','mcords3','num_chan',16,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0');  
param.max_tx = 30000; param.max_data_rate = 135; param.flight_hours = 7; param.sys_delay = 12.24e-6;
param.max_duty_cycle = 0.2;
param.create_IQ = false;
param.tg.altitude_guard = 750*12*2.54/100;
param.tg.staged_recording = true;
param.tg.Haltitude = 1250*12*2.54/100;
param.tg.Hice_thick = 3500;
param.fn = fullfile(base_dir,'noise_mode.xml');
param.prf = 12000;
param.presums = [3 3 29];
param.wfs(1).atten = 30;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
DDS_amp = final_DDS_amp;
DDS_amp = round(DDS_amp * param.max_tx / max(DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.1;
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 3e-6;
param.wfs(3).Tpd = 10e-6;
param.phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
[param.wfs(:).tx_mask] = deal([1 1 1 1 1 1 1 1]);
write_cresis_xml(param);

%% Image Mode
% All Ice <2250 m, 20000 +/-1500 ft AGL
param = struct('radar_name','mcords3','num_chan',16,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0');  
param.max_tx = 30000; param.max_data_rate = 300; param.flight_hours = 7; param.sys_delay = 12.24e-6;
param.max_duty_cycle = 0.2;
param.create_IQ = false;
param.tg.altitude_guard = 1500*12*2.54/100;
param.tg.staged_recording = false;
param.tg.Haltitude = 12000*12*2.54/100;
param.tg.Hice_thick = 2500;
param.fn = fullfile(base_dir,'image_mode_10us_3wf.xml');
param.prf = 10000;
param.presums = [13 11 13];
param.wfs(1).atten = 0;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
DDS_amp = final_DDS_amp .* Hwindow35 ./ Hwindow_orig;
DDS_amp = round(DDS_amp * param.max_tx / max(DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.1;
param.Tpd = 10e-6;
param.phase = final_DDS_phase;
param.tg.look_angle = [-30 0 30];
for wf = 1:length(param.tg.look_angle)
  nadir_vec = [0 0 1];
  beam_vec = [0 sind(param.tg.look_angle(wf)) cosd(param.tg.look_angle(wf))];
  nadir_delay = -nadir_vec * phase_centers;
  beam_delay = -beam_vec * phase_centers;
  nadir_delay = nadir_delay - nadir_delay(4);
  beam_delay = beam_delay - beam_delay(4);
  nadir_delay = nadir_delay / c;
  beam_delay = beam_delay / c - nadir_delay; % Only include delay relative to nadir
  beam_delay = [beam_delay - mean(beam_delay) 0];
  param.wfs(wf).delay = (final_DDS_time/1e9 - beam_delay)*1e9;
  % To do phase only correction:
  %param.wfs(wf).phase = final_DDS_phase + beam_delay * 2*pi*(f0+f1)/2;
  %param.wfs(wf).delay = final_DDS_time;
end
param.wfs(1).f0 = f0;
param.wfs(1).f1 = f1;
param.wfs(2).f0 = f0;
param.wfs(2).f1 = f1;
param.wfs(3).f0 = f0;
param.wfs(3).f1 = f1;
[param.wfs(:).tx_mask] = deal([0 0 0 0 0 0 0 0]);
write_cresis_xml(param);

%% Image Mode
% All Ice <1000 m, 1000 +/- 500 m AGL
param = struct('radar_name','mcords3','num_chan',16,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0');  
param.max_tx = 30000; param.max_data_rate = 300; param.flight_hours = 7; param.sys_delay = 12.24e-6;
param.max_duty_cycle = 0.2;
param.create_IQ = false;
param.tg.altitude_guard = 500;
param.tg.staged_recording = false;
param.tg.Haltitude = 1000;
param.tg.Hice_thick = 1250;
param.fn = fullfile(base_dir,'image_mode_3us_3wf.xml');
param.prf = 12000;
param.presums = [13 13 13];
param.wfs(1).atten = 10;
param.wfs(2).atten = 15;
param.wfs(3).atten = 10;
DDS_amp = final_DDS_amp .* [hanning(7).', 0] ./ Hwindow_orig;
DDS_amp = round(DDS_amp * param.max_tx / max(DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.1;
param.Tpd = 3e-6;
param.phase = final_DDS_phase;
param.tg.look_angle = [-30 0 30];
for wf = 1:length(param.tg.look_angle)
  nadir_vec = [0 0 1];
  beam_vec = [0 sind(param.tg.look_angle(wf)) cosd(param.tg.look_angle(wf))];
  nadir_delay = -nadir_vec * phase_centers;
  beam_delay = -beam_vec * phase_centers;
  nadir_delay = nadir_delay - nadir_delay(4);
  beam_delay = beam_delay - beam_delay(4);
  nadir_delay = nadir_delay / c;
  beam_delay = beam_delay / c - nadir_delay; % Only include delay relative to nadir
  beam_delay = [beam_delay - mean(beam_delay) 0];
  param.wfs(wf).delay = (final_DDS_time/1e9 - beam_delay)*1e9;
end
param.wfs(1).f0 = f0;
param.wfs(1).f1 = f1;
param.wfs(2).f0 = f0;
param.wfs(2).f1 = f1;
param.wfs(3).f0 = f0;
param.wfs(3).f1 = f1;
[param.wfs(:).tx_mask] = deal([0 0 0 0 0 0 0 0]);
write_cresis_xml(param);

%% Image Mode
% All Ice <1000 m, 500 +/- 500 m AGL
param = struct('radar_name','mcords3','num_chan',16,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0');  
param.max_tx = 30000; param.max_data_rate = 300; param.flight_hours = 7; param.sys_delay = 12.24e-6;
param.max_duty_cycle = 0.2;
param.create_IQ = false;
param.tg.altitude_guard = 500;
param.tg.staged_recording = false;
param.tg.Haltitude = 500;
param.tg.Hice_thick = 1250;
param.fn = fullfile(base_dir,'image_mode_1us_3wf.xml');
param.prf = 12000;
param.presums = [13 13 13];
param.wfs(1).atten = 15;
param.wfs(2).atten = 20;
param.wfs(3).atten = 15;
DDS_amp = final_DDS_amp .* [hanning(7).', 0] ./ Hwindow_orig;
DDS_amp = round(DDS_amp * param.max_tx / max(DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.2;
param.Tpd = 1e-6;
param.phase = final_DDS_phase;
param.tg.look_angle = [-30 0 30];
for wf = 1:length(param.tg.look_angle)
  nadir_vec = [0 0 1];
  beam_vec = [0 sind(param.tg.look_angle(wf)) cosd(param.tg.look_angle(wf))];
  nadir_delay = -nadir_vec * phase_centers;
  beam_delay = -beam_vec * phase_centers;
  nadir_delay = nadir_delay - nadir_delay(4);
  beam_delay = beam_delay - beam_delay(4);
  nadir_delay = nadir_delay / c;
  beam_delay = beam_delay / c - nadir_delay; % Only include delay relative to nadir
  beam_delay = [beam_delay - mean(beam_delay) 0];
  param.wfs(wf).delay = (final_DDS_time/1e9 - beam_delay)*1e9;
end
param.wfs(1).f0 = f0;
param.wfs(1).f1 = f1;
param.wfs(2).f0 = f0;
param.wfs(2).f1 = f1;
param.wfs(3).f0 = f0;
param.wfs(3).f1 = f1;
[param.wfs(:).tx_mask] = deal([0 0 0 0 0 0 0 0]);
write_cresis_xml(param);

%% Image Mode
% All Ice <1500 m, 1000 +/- 500 m AGL
param = struct('radar_name','mcords3','num_chan',16,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0');  
param.max_tx = 30000; param.max_data_rate = 300; param.flight_hours = 7; param.sys_delay = 12.24e-6;
param.max_duty_cycle = 0.2;
param.create_IQ = false;
param.tg.altitude_guard = 500;
param.tg.staged_recording = false;
param.tg.Haltitude = 1000;
param.tg.Hice_thick = 1750;
param.fn = fullfile(base_dir,'image_mode_sarqardliupsermia.xml');
param.prf = 12000;
param.presums = [13 13 13];
param.wfs(1).atten = 0;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
DDS_amp = final_DDS_amp .* [hanning(7).', 0] ./ Hwindow_orig;
DDS_amp = round(DDS_amp * param.max_tx / max(DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.1;
param.Tpd = 3e-6;
param.phase = final_DDS_phase;
param.tg.look_angle = [-30 0 30];
for wf = 1:length(param.tg.look_angle)
  nadir_vec = [0 0 1];
  beam_vec = [0 sind(param.tg.look_angle(wf)) cosd(param.tg.look_angle(wf))];
  nadir_delay = -nadir_vec * phase_centers;
  beam_delay = -beam_vec * phase_centers;
  nadir_delay = nadir_delay - nadir_delay(4);
  beam_delay = beam_delay - beam_delay(4);
  nadir_delay = nadir_delay / c;
  beam_delay = beam_delay / c - nadir_delay; % Only include delay relative to nadir
  beam_delay = [beam_delay - mean(beam_delay) 0];
  param.wfs(wf).delay = (final_DDS_time/1e9 - beam_delay)*1e9;
end
param.wfs(1).f0 = f0;
param.wfs(1).f1 = f1;
param.wfs(2).f0 = f0;
param.wfs(2).f1 = f1;
param.wfs(3).f0 = f0;
param.wfs(3).f1 = f1;
[param.wfs(:).tx_mask] = deal([0 0 0 0 0 0 0 0]);
write_cresis_xml(param);

%% Image Mode
% All Ice <2250 m, 20000 +/-1500 ft AGL
param = struct('radar_name','mcords3','num_chan',16,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0');  
param.max_tx = 30000; param.max_data_rate = 350; param.flight_hours = 7; param.sys_delay = 12.24e-6;
param.max_duty_cycle = 0.2;
param.create_IQ = false;
param.tg.altitude_guard = 1500*12*2.54/100;
param.tg.staged_recording = false;
param.tg.Haltitude = 12000*12*2.54/100;
param.tg.Hice_thick = 2500;
param.fn = fullfile(base_dir,'image_mode_10us_5wf.xml');
param.prf = 10000;
param.presums = [7 5 5 5 5];
param.wfs(1).atten = 0;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
param.wfs(4).atten = 0;
param.wfs(5).atten = 0;
DDS_amp = final_DDS_amp .* [0.5 1 1 1 1 1 0.5 0] ./ Hwindow_orig;
DDS_amp = round(DDS_amp * param.max_tx / max(DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.1; 
param.Tpd = 10e-6;
param.phase = final_DDS_phase;
% Set time delays for beam forming
param.tg.look_angle = [-34 -17 0 17 34];
for wf = 1:length(param.tg.look_angle)
  nadir_vec = [0 0 1];
  beam_vec = [0 sind(param.tg.look_angle(wf)) cosd(param.tg.look_angle(wf))];
  nadir_delay = -nadir_vec * phase_centers;
  beam_delay = -beam_vec * phase_centers;
  nadir_delay = nadir_delay - nadir_delay(4);
  beam_delay = beam_delay - beam_delay(4);
  nadir_delay = nadir_delay / c;
  beam_delay = beam_delay / c - nadir_delay; % Only include delay relative to nadir
  beam_delay = [beam_delay - mean(beam_delay) 0];
  param.wfs(wf).delay = (final_DDS_time/1e9 - beam_delay)*1e9;
end
param.wfs(1).f0 = f0;
param.wfs(1).f1 = f1;
param.wfs(2).f0 = f0;
param.wfs(2).f1 = f1;
param.wfs(3).f0 = f0;
param.wfs(3).f1 = f1;
param.wfs(4).f0 = f0;
param.wfs(4).f1 = f1;
param.wfs(5).f0 = f0;
param.wfs(5).f1 = f1;
[param.wfs(:).tx_mask] = deal([0 0 0 0 0 0 0 0]);
write_cresis_xml(param);

param.presums = [15 15 15 15 15];
param.wfs(1).atten = 15;
param.wfs(2).atten = 15;
param.wfs(3).atten = 15;
param.wfs(4).atten = 15;
param.wfs(5).atten = 15;
param.fn = fullfile(base_dir,'image_mode_10us_5wf_tukey_roll_test.xml');
write_cresis_xml(param);

param.presums = [15 15 15 15 15];
param.wfs(1).atten = 10;
param.wfs(2).atten = 10;
param.wfs(3).atten = 10;
param.wfs(4).atten = 10;
param.wfs(5).atten = 10;
DDS_amp = final_DDS_amp .* [hanning(7)' 0] ./ Hwindow_orig;
param.fn = fullfile(base_dir,'image_mode_10us_5wf_hanning_roll_test.xml');
write_cresis_xml(param);


%% Equalization High Altitude (Using Ocean)
% 20000 ft +/- 2500 ft
% For lower altitude, increase attenuation
% Use these settings over ocean or sea ice for fast-time equalization,
% transmit equalization, and receiver equalization
for Tpd = [1e-6 3e-6 10e-6]
  param = struct('radar_name','mcords3','num_chan',16,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0');  
  param.max_tx = 30000; param.max_data_rate = 135; param.flight_hours = 7; param.sys_delay = 12.24e-6;
  param.max_duty_cycle = 0.2;
  param.create_IQ = false;
  param.tg.altitude_guard = 6000*12*2.54/100;
  param.tg.staged_recording = false;
  param.tg.Haltitude = 19000*12*2.54/100;
  param.tg.Hice_thick = 0;
  param.fn = fullfile(base_dir,sprintf('equalization_%.0fus.xml', Tpd*1e6));
  param.prf = 10000;
  param.presums = [15 15 15 15 15 15 15 15];
  [param.wfs(1:8).atten] = deal(5);
  param.wfs(8).atten = deal(15);
  param.tx_weights = final_DDS_amp;
  if Tpd > 5e-6
    param.tukey = 0.1;
  else
    param.tukey = 0.2;
  end
  param.Tpd = Tpd;
  param.phase = final_DDS_phase;
  param.delay = final_DDS_time;
  param.f0 = f0;
  param.f1 = f1;
  for wf=1:7
    param.wfs(wf).tx_mask = ones(1,8);
    param.wfs(wf).tx_mask(9-wf) = 0;
  end
  param.wfs(8).tx_mask = zeros(1,8);
  write_cresis_xml(param);
end
