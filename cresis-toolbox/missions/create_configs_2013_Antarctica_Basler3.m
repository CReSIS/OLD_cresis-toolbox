% script create_configs_2013_Antarctica_Basler3
%
% Creates NI radar depth sounder settings

base_dir = '/basler/scratch2/';
base_dir = 'c:/tmp/';
final_DDS_phase = [-43	-21	60.3	-45.7	-27.4	-62.4	-25.6	-65.6];
final_DDS_amp = [22118	33496	36455	29475	24741	26106	50000	12948];
final_DDS_time = [-0.51	-0.27	0.68	-0.61	-0.54	-0.69	-0.53	-0.87];
f0 = 200e6;
f1 = 450e6;

Hwindow_orig = [0.5 1 1 1 1 1 1 0.5]; % Window created by final_DDS_amp weights
Hwindow35 = [0.1915	0.4636	0.7843	1	1	0.7843	0.4636	0.1915]; % chebwin(8,35)
Hwindow55 = [0.0795	0.324	0.7014	1	1	0.7014	0.324	0.0795];

%% Survey Mode
% <1500 m
param = struct('radar_name','mcords4','num_chan',8,'aux_dac',[192 192 192 192 192 192 192 192],'version','12.0.1f3');
param.max_tx = 50000; param.max_data_rate = 600; param.flight_hours = 7;
param.max_duty_cycle = 0.16;
param.create_IQ = true;
param.tg.altitude_guard = 750*12*2.54/100;
param.tg.staged_recording = true;
param.tg.Haltitude = 1000*12*2.54/100;
param.tg.Hice_thick = 1500;
param.fn = fullfile(base_dir,'survey_mode_3us_2wf_1500mthick.xml');
param.prf = 12000;
param.presums = [3 33];
param.wfs(1).atten = 25;
param.wfs(2).atten = 0;
DDS_amp = final_DDS_amp .* Hwindow35 ./ Hwindow_orig;
DDS_amp = round(DDS_amp * param.max_tx / max(DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.03;
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 3e-6;
param.phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
[param.wfs(:).tx_mask] = deal(0);
write_cresis_xml(param);
param.presums = [3 17];
param.fn = fullfile(base_dir,'survey_mode_3us_2wf_1500mthick_higheprf.xml');
write_cresis_xml(param);
param.presums = [3 33];
DDS_amp = param.max_tx * ones(size(final_DDS_amp));
param.tx_weights = DDS_amp;
param.fn = fullfile(base_dir,'survey_mode_3us_2wf_1500mthick_maxpower.xml');
write_cresis_xml(param);
param.presums = [3 17];
param.fn = fullfile(base_dir,'survey_mode_3us_2wf_1500mthick_maxpower_higheprf.xml');
write_cresis_xml(param);

%% Survey Mode
% >800 m and <2000 m
param = struct('radar_name','mcords4','num_chan',8,'aux_dac',[192 192 192 192 192 192 192 192],'version','12.0.1f3'); 
param.max_tx = 50000; param.max_data_rate = 600; param.flight_hours = 7;
param.max_duty_cycle = 0.16;
param.create_IQ = true;
param.tg.altitude_guard = 750*12*2.54/100;
param.tg.staged_recording = true;
param.tg.Haltitude = 1000*12*2.54/100;
param.tg.Hice_thick = 2000;
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_2000mthick.xml');
param.prf = 12000;
param.presums = [3 3 33];
param.wfs(1).atten = 25;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
DDS_amp = final_DDS_amp .* Hwindow35 ./ Hwindow_orig;
DDS_amp = round(DDS_amp * param.max_tx / max(DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.03; 
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 3e-6;
param.wfs(3).Tpd = 10e-6;
param.phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
[param.wfs(:).tx_mask] = deal(0);
write_cresis_xml(param);
param.presums = [3 3 17];
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_2000mthick_higheprf.xml');
write_cresis_xml(param);
param.presums = [3 3 33];
DDS_amp = param.max_tx * ones(size(final_DDS_amp));
param.tx_weights = DDS_amp;
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_2000mthick_maxpower.xml');
write_cresis_xml(param);

%% Survey Mode
% >1500 m and <3000 m
param = struct('radar_name','mcords4','num_chan',8,'aux_dac',[192 192 192 192 192 192 192 192],'version','12.0.1f3'); 
param.max_tx = 50000; param.max_data_rate = 600; param.flight_hours = 7;
param.max_duty_cycle = 0.16;
param.create_IQ = true;
param.tg.altitude_guard = 750*12*2.54/100;
param.tg.staged_recording = true;
param.tg.Haltitude = 1000*12*2.54/100;
param.tg.Hice_thick = 3000;
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_3000mthick.xml');
param.prf = 12000;
param.presums = [3 3 33];
param.wfs(1).atten = 25;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
DDS_amp = final_DDS_amp .* Hwindow35 ./ Hwindow_orig;
DDS_amp = round(DDS_amp * param.max_tx / max(DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.03; 
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 3e-6;
param.wfs(3).Tpd = 10e-6;
param.phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
[param.wfs(:).tx_mask] = deal(0);
write_cresis_xml(param);
param.presums = [3 3 17];
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_3000mthick_higheprf.xml');
write_cresis_xml(param);
param.presums = [3 3 33];
DDS_amp = param.max_tx * ones(size(final_DDS_amp));
param.tx_weights = DDS_amp;
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_3000mthick_maxpower.xml');
write_cresis_xml(param);

%% Survey Mode
% >3000 m
param = struct('radar_name','mcords4','num_chan',8,'aux_dac',[192 192 192 192 192 192 192 192],'version','12.0.1f3'); 
param.max_tx = 50000; param.max_data_rate = 600; param.flight_hours = 7;
param.max_duty_cycle = 0.16;
param.create_IQ = true;
param.tg.altitude_guard = 750*12*2.54/100;
param.tg.staged_recording = true;
param.tg.Haltitude = 1000*12*2.54/100;
param.tg.Hice_thick = 4500;
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_4500mthick.xml');
param.prf = 12000;
param.presums = [3 3 33];
param.wfs(1).atten = 25;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
DDS_amp = final_DDS_amp .* Hwindow35 ./ Hwindow_orig;
DDS_amp = round(DDS_amp * param.max_tx / max(DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.03; 
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 3e-6;
param.wfs(3).Tpd = 10e-6;
param.phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
[param.wfs(:).tx_mask] = deal(0);
write_cresis_xml(param);
param.presums = [3 3 17];
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_4500mthick_higheprf.xml');
write_cresis_xml(param);
param.presums = [3 3 33];
DDS_amp = param.max_tx * ones(size(final_DDS_amp));
param.tx_weights = DDS_amp;
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_4500mthick_maxpower.xml');
write_cresis_xml(param);
param.presums = [3 3 17];
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_4500mthick_maxpower_higheprf.xml');
write_cresis_xml(param);

%% Survey Mode (Noise Measurement)
% NA
param = struct('radar_name','mcords4','num_chan',8,'aux_dac',[192 192 192 192 192 192 192 192],'version','12.0.1f3'); 
param.max_tx = 50000; param.max_data_rate = 600; param.flight_hours = 7;
param.max_duty_cycle = 0.16;
param.create_IQ = true;
param.tg.altitude_guard = 750*12*2.54/100;
param.tg.staged_recording = true;
param.tg.Haltitude = 1000*12*2.54/100;
param.tg.Hice_thick = 1500;
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_noise.xml');
param.prf = 12000;
param.presums = [3 3 33];
param.wfs(1).atten = 25;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
DDS_amp = 0 * ones(size(final_DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.03; 
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 3e-6;
param.wfs(3).Tpd = 10e-6;
param.phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
[param.wfs(:).tx_mask] = deal([1 1 1 1 1 1 1 1]);
write_cresis_xml(param);

%% Survey Mode MAX POWER
% DEEPEST ICE
param = struct('radar_name','mcords4','num_chan',8,'aux_dac',[192 192 192 192 192 192 192 192],'version','12.0.1f3'); 
param.max_tx = 50000; param.max_data_rate = 600; param.flight_hours = 7;
param.max_duty_cycle = 0.16;
param.create_IQ = true;
param.tg.altitude_guard = 750*12*2.54/100;
param.tg.staged_recording = true;
param.tg.Haltitude = 1000*12*2.54/100;
param.tg.Hice_thick = 4500;
param.fn = fullfile(base_dir,'survey_mode_10us_2wf_narrowBW_maxpower.xml');
param.prf = 12000;
param.presums = [3 33];
param.wfs(1).atten = 10;
param.wfs(2).atten = 0;
DDS_amp = param.max_tx * ones(size(final_DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.03;
param.wfs(1).Tpd = 3e-6;
param.wfs(2).Tpd = 10e-6;
param.phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = 200e6;
param.f1 = 230e6;
[param.wfs(:).tx_mask] = deal(0);
write_cresis_xml(param);
param.presums = [3 17];
param.fn = fullfile(base_dir,'survey_mode_10us_2wf_narrowBW_maxpower_higheprf.xml');
write_cresis_xml(param);

%% Glacier Mode
% ? m
param = struct('radar_name','mcords4','num_chan',8,'aux_dac',[192 192 192 192 192 192 192 192],'version','12.0.1f3'); 
param.max_tx = 50000; param.max_data_rate = 600; param.flight_hours = 7;
param.max_duty_cycle = 0.16;
param.create_IQ = true;
param.tg.altitude_guard = 750*12*2.54/100;
param.tg.staged_recording = true;
param.tg.Haltitude = 1000*12*2.54/100;
param.tg.Hice_thick = 3000;
param.fn = fullfile(base_dir,'glacier_mode_2wf.xml');
param.prf = 12000;
param.presums = [3 17];
param.wfs(1).atten = 25;
param.wfs(2).atten = 0;
DDS_amp = final_DDS_amp .* Hwindow35 ./ Hwindow_orig;
DDS_amp = round(DDS_amp * param.max_tx / max(DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.03;
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 10e-6;
param.phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = 200e6;
param.f1 = 450e6;
[param.wfs(:).tx_mask] = deal(0);
write_cresis_xml(param);

%% Image Mode
% All Ice <1000 m, 6000 +/-1000 ft AGL
param = struct('radar_name','mcords4','num_chan',8,'aux_dac',[192 192 192 192 192 192 192 192],'version','12.0.1f3'); 
param.max_tx = 50000; param.max_data_rate = 350; param.flight_hours = 7;
param.max_duty_cycle = 0.16;
param.create_IQ = true;
param.tg.altitude_guard = 750*12*2.54/100;
param.tg.staged_recording = false;
param.tg.Haltitude = 6000*12*2.54/100;
param.tg.Hice_thick = 1000 + 150;
param.fn = fullfile(base_dir,'image_mode_3us_3wf.xml');
param.prf = 12000;
param.presums = [9 9 9];
param.wfs(1).atten = 0;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
DDS_amp = final_DDS_amp .* Hwindow35 ./ Hwindow_orig;
DDS_amp = round(DDS_amp * param.max_tx / max(DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.03; 
param.Tpd = 3e-6;
param.phase = final_DDS_phase;
param.wfs(1).delay = final_DDS_time;
param.wfs(2).delay = final_DDS_time - 0.55*(0:7);
param.wfs(3).delay = final_DDS_time + 0.55*(0:7);
param.f0 = f0;
param.f1 = f1;
[param.wfs(:).tx_mask] = deal(0);
write_cresis_xml(param);

%% Image Mode
% All Ice 1000-2500 m, 6000 +/-1000 ft AGL
param = struct('radar_name','mcords4','num_chan',8,'aux_dac',[192 192 192 192 192 192 192 192],'version','12.0.1f3'); 
param.max_tx = 50000; param.max_data_rate = 600; param.flight_hours = 7;
param.max_duty_cycle = 0.16;
param.create_IQ = true;
param.tg.altitude_guard = 1000*12*2.54/100;
param.tg.staged_recording = false;
param.tg.Haltitude = 6000*12*2.54/100;
param.tg.Hice_thick = 2500;
param.fn = fullfile(base_dir,'image_mode_10us_3wf.xml');
param.prf = 12000;
param.presums = [13 13 13];
param.wfs(1).atten = 0;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
DDS_amp = final_DDS_amp .* Hwindow35 ./ Hwindow_orig;
DDS_amp = round(DDS_amp * param.max_tx / max(DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.03; 
param.Tpd = 10e-6;
param.phase = final_DDS_phase;
param.wfs(1).delay = final_DDS_time;
param.wfs(2).delay = final_DDS_time - 0.55*(0:7);
param.wfs(3).delay = final_DDS_time + 0.55*(0:7);
param.f0 = f0;
param.f1 = f1;
[param.wfs(:).tx_mask] = deal(0);
write_cresis_xml(param);

%% Image (MIMO) Mode
% Ice >1000 m, 6000 +/-1000 ft AGL
param = struct('radar_name','mcords4','num_chan',8,'aux_dac',[192 192 192 192 192 192 192 192],'version','12.0.1f3'); 
param.max_tx = 50000; param.max_data_rate = 600; param.flight_hours = 7;
param.max_duty_cycle = 0.16;
param.create_IQ = true;
param.tg.altitude_guard = 1000*12*2.54/100;
param.tg.staged_recording = false;
param.tg.Haltitude = 6000*12*2.54/100;
param.tg.Hice_thick = 1000 + 750;
param.fn = fullfile(base_dir,'mimo_mode_3us_3wf.xml');
param.prf = 12000;
param.presums = [11 11];
param.wfs(1).atten = 0;
param.wfs(2).atten = 0;
param.tx_weights = [50000 50000 0 0 0 0 50000 50000];
param.tukey = 0.03; 
param.Tpd = 3e-6;
param.phase = final_DDS_phase;
param.delay = final_DDS_time;
param.wfs(1).tx_mask = [0 0 1 1 1 1 1 1];
param.wfs(2).tx_mask = [1 1 1 1 1 1 0 0];
param.wfs(1).f0 = f0;
param.wfs(1).f1 = f1;
param.wfs(2).f0 = 1e9-f1;
param.wfs(2).f1 = 1e9-f0;
write_cresis_xml(param);

%% Image (MIMO) Mode
% Ice >1000 m, 6000 +/-1000 ft AGL
param = struct('radar_name','mcords4','num_chan',8,'aux_dac',[192 192 192 192 192 192 192 192],'version','12.0.1f3'); 
param.max_tx = 50000; param.max_data_rate = 600; param.flight_hours = 7;
param.max_duty_cycle = 0.16;
param.create_IQ = true;
param.tg.altitude_guard = 1000*12*2.54/100;
param.tg.staged_recording = false;
param.tg.Haltitude = 6000*12*2.54/100;
param.tg.Hice_thick = 2500;
param.fn = fullfile(base_dir,'mimo_mode_10us_3wf.xml');
param.prf = 12000;
param.presums = [15 15];
param.wfs(1).atten = 0;
param.wfs(2).atten = 0;
param.tx_weights = [50000 50000 0 0 0 0 50000 50000];
param.tukey = 0.03; 
param.Tpd = 10e-6;
param.phase = final_DDS_phase;
param.delay = final_DDS_time;
param.wfs(1).tx_mask = [0 0 1 1 1 1 1 1];
param.wfs(2).tx_mask = [1 1 1 1 1 1 0 0];
param.wfs(1).f0 = f0;
param.wfs(1).f1 = f1;
param.wfs(2).f0 = f1;
param.wfs(2).f1 = f0;
write_cresis_xml(param);

%% Transmit Equalization Low Altitude (Using Flat Bed)
% 1000 ft +/- 500 ft
param = struct('radar_name','mcords4','num_chan',8,'aux_dac',[192 192 192 192 192 192 192 192],'version','12.0.1f3'); 
param.max_tx = 50000; param.max_data_rate = 600; param.flight_hours = 7;
param.max_duty_cycle = 0.16;
param.create_IQ = true;
param.tg.altitude_guard = 750*12*2.54/100;
param.tg.staged_recording = false;
param.tg.Haltitude = 1000*12*2.54/100;
param.tg.Hice_thick = 1000;
param.fn = fullfile(base_dir,'tx_equal_3us_16wf_lowAGL.xml');
param.prf = 12000;
param.presums = [15 15 15 15 15 15 15 15];
[param.wfs(1:8).atten] = deal(0);
param.tx_weights = final_DDS_amp;
param.tukey = 0.03; 
param.Tpd = 3e-6;
param.phase = final_DDS_phase;
param.time = final_DDS_time;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
for wf=1:8
  param.wfs(wf).tx_mask = ones(1,8);
  param.wfs(wf).tx_mask(wf) = 0;
end
write_cresis_xml(param);

%% Transmit Equalization High Altitude
% 8000 +/-2000 ft AGL
param = struct('radar_name','mcords4','num_chan',8,'aux_dac',[192 192 192 192 192 192 192 192],'version','12.0.1f3'); 
param.max_tx = 50000; param.max_data_rate = 600; param.flight_hours = 7;
param.max_duty_cycle = 0.16;
param.create_IQ = true;
param.tg.altitude_guard = 3000*12*2.54/100;
param.tg.staged_recording = false;
param.tg.Haltitude = 8000*12*2.54/100;
param.tg.Hice_thick = 0;
param.fn = fullfile(base_dir,'tx_equal_1us_16wf_highAGL.xml');
param.prf = 12000;
param.presums = [15 15 15 15 15 15 15 15];
[param.wfs(1:8).atten] = deal(15);
param.tx_weights = final_DDS_amp;
param.tukey = 0.03; 
param.Tpd = 1e-6;
param.phase = final_DDS_phase;
param.time = final_DDS_time;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
for wf=1:8
  param.wfs(wf).tx_mask = ones(1,8);
  param.wfs(wf).tx_mask(wf) = 0;
end
write_cresis_xml(param);

%% Transmit Equalization High Altitude
% 8000 +/-2000 ft AGL
param = struct('radar_name','mcords4','num_chan',8,'aux_dac',[192 192 192 192 192 192 192 192],'version','12.0.1f3'); 
param.max_tx = 50000; param.max_data_rate = 600; param.flight_hours = 7;
param.max_duty_cycle = 0.16;
param.create_IQ = true;
param.tg.altitude_guard = 3000*12*2.54/100;
param.tg.staged_recording = false;
param.tg.Haltitude = 8000*12*2.54/100;
param.tg.Hice_thick = 0;
param.fn = fullfile(base_dir,'tx_equal_3us_16wf_highAGL.xml');
param.prf = 12000;
param.presums = [15 15 15 15 15 15 15 15];
[param.wfs(1:8).atten] = deal(15);
param.tx_weights = final_DDS_amp;
param.tukey = 0.03; 
param.Tpd = 3e-6;
param.phase = final_DDS_phase;
param.time = final_DDS_time;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
for wf=1:8
  param.wfs(wf).tx_mask = ones(1,8);
  param.wfs(wf).tx_mask(wf) = 0;
end
write_cresis_xml(param);

%% Transmit Equalization High Altitude
% 8000 +/-2000 ft AGL
param = struct('radar_name','mcords4','num_chan',8,'aux_dac',[192 192 192 192 192 192 192 192],'version','12.0.1f3'); 
param.max_tx = 50000; param.max_data_rate = 600; param.flight_hours = 7;
param.max_duty_cycle = 0.16;
param.create_IQ = true;
param.tg.altitude_guard = 3000*12*2.54/100;
param.tg.staged_recording = false;
param.tg.Haltitude = 8000*12*2.54/100;
param.tg.Hice_thick = 0;
param.fn = fullfile(base_dir,'tx_equal_10us_16wf_highAGL.xml');
param.prf = 12000;
param.presums = [15 15 15 15 15 15 15 15];
[param.wfs(1:8).atten] = deal(15);
param.tx_weights = final_DDS_amp;
param.tukey = 0.03; 
param.Tpd = 10e-6;
param.phase = final_DDS_phase;
param.time = final_DDS_time;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
for wf=1:8
  param.wfs(wf).tx_mask = ones(1,8);
  param.wfs(wf).tx_mask(wf) = 0;
end
write_cresis_xml(param);

%% Transmit Beam Pattern High Altitude
% 8000 +/-2000 ft AGL
param = struct('radar_name','mcords4','num_chan',8,'aux_dac',[192 192 192 192 192 192 192 192],'version','12.0.1f3'); 
param.max_tx = 50000; param.max_data_rate = 600; param.flight_hours = 7;
param.max_duty_cycle = 0.16;
param.create_IQ = true;
param.tg.altitude_guard = 3000*12*2.54/100;
param.tg.staged_recording = false;
param.tg.Haltitude = 8000*12*2.54/100;
param.tg.Hice_thick = 1000;
param.fn = fullfile(base_dir,'tx_beam_pattern_10us_6wf_highAGL.xml');
param.prf = 12000;
param.presums = [15 15 15 15 15 15 15 15];
[param.wfs(1:3).atten] = deal(30);
DDS_amp = final_DDS_amp .* Hwindow35 ./ Hwindow_orig;
DDS_amp = round(DDS_amp * param.max_tx / max(DDS_amp));
param.tx_weights = DDS_amp;
param.tukey = 0.03; 
param.Tpd = 10e-6;
param.phase = final_DDS_phase;
param.wfs(1).delay = final_DDS_time;
param.wfs(2).delay = final_DDS_time - 0.55*(0:7);
param.wfs(3).delay = final_DDS_time + 0.55*(0:7);
param.f0 = f0;
param.f1 = f1;
[param.wfs(:).tx_mask] = deal(0);
write_cresis_xml(param);
