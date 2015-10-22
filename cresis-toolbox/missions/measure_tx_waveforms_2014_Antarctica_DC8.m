param = [];
param.num_ave = 256;
param.num_pnts = 1000000;
param.fs = 5.0e9;
param.trigger = 'CH3';
param.channels = {'CH1','CH2','CH3'};
% 500 MHz, 2.5 GSPS DSO (see Utility->IO->USB Computer->USBTMC for address)
%   param.visa_address = 'USB0::0x0699::0x0401::C000466::0::INSTR';
% 1 GHz, 5 GSPS DSO (see Utility->IO->USB Computer->USBTMC for address)
param.visa_address = 'USB0::0x0699::0x0401::C021536::0::INSTR';
param.time_out_sec = 60;
[instrument,time,data] = tektronix_dso_capture(param);

figure(1); clf;
plot(time,data(:,1),'y');
hold on;
plot(time,data(:,3),'m');
plot(time,data(:,2),'b');
hold off;
grid on;

% Spare EPRI cable was used to route EPRI signal from DDS to DSO chan 3
% Another random cable (6? ft) used to route TTL0 to DSO chan 1 (for DDS output measurements)
% Another random cable (6? ft) used to route TTL2 to DSO chan 1 (for power Amp output measurements)
% TTL0-TTL7 are all the same timing

% param.counts = 15000;
% notes = sprintf('%s: Output of TxAmp with %d counts --> Type-N cable --> 30 dB high power atten --> TypeN/SMA adapter --> 20 dBm atten --> tan GPS-1 SMA cable --> CH2-DSO, CH1-TTL0, CH3-EPRI (all other DDS outputs 50ohm)', datestr(now), param.counts);
% save('txamp_chan1_10000ns.mat','time','data','notes','param');
% save('txamp_chan1_3000ns.mat','time','data','notes','param');
% save('txamp_chan1_1000ns.mat','time','data','notes','param');
% save('txamp_chan1_30ns.mat','time','data','notes','param');
% 
% param.counts = 15000;
% notes = sprintf('%s: Output of TxAmp with %d counts --> Type-N cable --> 30 dB high power atten --> TypeN/SMA adapter --> 20 dBm atten --> tan GPS-1 SMA cable --> CH2-DSO, CH1-TTL0, CH3-EPRI (all other DDS outputs 50ohm)', datestr(now), param.counts);
% save('txamp_chan2_10000ns.mat','time','data','notes','param');
% save('txamp_chan2_3000ns.mat','time','data','notes','param');
% save('txamp_chan2_1000ns.mat','time','data','notes','param');
% save('txamp_chan2_30ns.mat','time','data','notes','param');
% 
% param.counts = 17000;
% notes = sprintf('%s: Output of TxAmp with %d counts --> Type-N cable --> 30 dB high power atten --> TypeN/SMA adapter --> 20 dBm atten --> tan GPS-1 SMA cable --> CH2-DSO, CH1-TTL0, CH3-EPRI (all other DDS outputs 50ohm)', datestr(now), param.counts);
% save('txamp_chan3_10000ns.mat','time','data','notes','param');
% save('txamp_chan3_3000ns.mat','time','data','notes','param');
% save('txamp_chan3_1000ns.mat','time','data','notes','param');
% save('txamp_chan3_30ns.mat','time','data','notes','param');
% 
% param.counts = 15920;
% notes = sprintf('%s: Output of TxAmp with %d counts --> Type-N cable --> 30 dB high power atten --> TypeN/SMA adapter --> 20 dBm atten --> tan GPS-1 SMA cable --> CH2-DSO, CH1-TTL0, CH3-EPRI (all other DDS outputs 50ohm)', datestr(now), param.counts);
% save('txamp_chan4_10000ns.mat','time','data','notes','param');
% save('txamp_chan4_3000ns.mat','time','data','notes','param');
% save('txamp_chan4_1000ns.mat','time','data','notes','param');
% save('txamp_chan4_30ns.mat','time','data','notes','param');
% 
% param.counts = 15000;
% notes = sprintf('%s: Output of TxAmp with %d counts --> Type-N cable --> 30 dB high power atten --> TypeN/SMA adapter --> 20 dBm atten --> tan GPS-1 SMA cable --> CH2-DSO, CH1-TTL0, CH3-EPRI (all other DDS outputs 50ohm)', datestr(now), param.counts);
% save('txamp_chan5_10000ns.mat','time','data','notes','param');
% save('txamp_chan5_3000ns.mat','time','data','notes','param');
% save('txamp_chan5_1000ns.mat','time','data','notes','param');
% save('txamp_chan5_30ns.mat','time','data','notes','param');
% 
% param.counts = 14150;
% notes = sprintf('%s: Output of TxAmp with %d counts --> Type-N cable --> 30 dB high power atten --> TypeN/SMA adapter --> 20 dBm atten --> tan GPS-1 SMA cable --> CH2-DSO, CH1-TTL0, CH3-EPRI (all other DDS outputs 50ohm)', datestr(now), param.counts);
% save('txamp_chan6_10000ns.mat','time','data','notes','param');
% save('txamp_chan6_3000ns.mat','time','data','notes','param');
% save('txamp_chan6_1000ns.mat','time','data','notes','param');
% save('txamp_chan6_30ns.mat','time','data','notes','param');


% Calibration 100-300 MHz
% 4 port NA N5222A, port 1 and port 2
% 
% Measurement 1
% Port 1 --> 36" UFA147A --> REFERENCE PLANE --> female SMA-female Type-N (added) --> Pasternack RG214 5? ft --> 30 dB 200 W atten --> 20? tan GPS-1 SMA cable (missing SMA female-BNC male adapter) --> REFERENCE PLANE --> fem-fem SMA adapter --> 36" UFA147A --> Port 2
% 
% Measurement 2
% Port 1 --> 36" UFA147A --> REFERENCE PLANE --> female SMA-female SMA (added) --> 20? tan GPS-1 SMA cable (missing SMA female-BNC male adapter) --> REFERENCE PLANE --> fem-fem SMA adapter --> 36" UFA147A --> Port 2
% 
% Measurement 3-8, RADAR = Antenna input on box to ADC cable output (right before ADC)
% Port 1 --> 36" UFA147A --> REFERENCE PLANE --> female SMA-female Type-N (added) --> Pasternack RG214 5? ft --> 30 dB 200 W atten --> 20? tan GPS-1 SMA cable (missing SMA female-BNC male adapter) --> SMA fem-TypeN male (added) -> RADAR --> REFERENCE PLANE --> fem-fem SMA adapter --> 36" UFA147A --> Port 2





% notes = sprintf('%s: Output of DDS with 17000 counts --> BNC/SMA adapter --> tan GPS-1 SMA cable --> CH2-DSO, CH1-TTL0, CH3-EPRI (all other DDS outputs 50ohm)', datestr(now));
% save('chan1_10000ns.mat','time','data','notes','param');
% save('chan1_3000ns.mat','time','data','notes','param');
% save('chan1_1000ns.mat','time','data','notes','param');
% save('chan1_30ns.mat','time','data','notes','param');
% 
% notes = sprintf('%s: Output of DDS with 17000 counts --> BNC/SMA adapter --> tan GPS-1 SMA cable --> CH2-DSO, CH1-TTL0, CH3-EPRI (all other DDS outputs open)', datestr(now));
% save('chan2_10000ns.mat','time','data','notes','param');
% save('chan2_3000ns.mat','time','data','notes','param');
% save('chan2_1000ns.mat','time','data','notes','param');
% save('chan2_30ns.mat','time','data','notes','param');
% 
% notes = sprintf('%s: Output of DDS with 17000 counts --> BNC/SMA adapter --> tan GPS-1 SMA cable --> CH2-DSO, CH1-TTL0, CH3-EPRI (all other DDS outputs 50ohm)', datestr(now));
% save('chan3_10000ns.mat','time','data','notes','param');
% save('chan3_3000ns.mat','time','data','notes','param');
% save('chan3_1000ns.mat','time','data','notes','param');
% save('chan3_30ns.mat','time','data','notes','param');
% 
% notes = sprintf('%s: Output of DDS with 17000 counts --> BNC/SMA adapter --> tan GPS-1 SMA cable --> CH2-DSO, CH1-TTL0, CH3-EPRI (all other DDS outputs 50ohm)', datestr(now));
% save('chan4_10000ns.mat','time','data','notes','param');
% save('chan4_3000ns.mat','time','data','notes','param');
% save('chan4_1000ns.mat','time','data','notes','param');
% save('chan4_30ns.mat','time','data','notes','param');
% 
% notes = sprintf('%s: Output of DDS with 17000 counts --> BNC/SMA adapter --> tan GPS-1 SMA cable --> CH2-DSO, CH1-TTL0, CH3-EPRI (all other DDS outputs 50ohm)', datestr(now));
% save('chan5_10000ns.mat','time','data','notes','param');
% save('chan5_3000ns.mat','time','data','notes','param');
% save('chan5_1000ns.mat','time','data','notes','param');
% save('chan5_30ns.mat','time','data','notes','param');
% 
% notes = sprintf('%s: Output of DDS with 17000 counts --> BNC/SMA adapter --> tan GPS-1 SMA cable --> CH2-DSO, CH1-TTL0, CH3-EPRI (all other DDS outputs 50ohm)', datestr(now));
% save('chan6_10000ns.mat','time','data','notes','param');
% save('chan6_3000ns.mat','time','data','notes','param');
% save('chan6_1000ns.mat','time','data','notes','param');
% save('chan6_30ns.mat','time','data','notes','param');










