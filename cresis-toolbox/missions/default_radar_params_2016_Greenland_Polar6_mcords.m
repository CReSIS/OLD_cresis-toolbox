% script default_radar_params_2016_Greenland_Polar6_mcords

default.radar.fs = 1600e6;
default.radar.adc_bits = 12;
default.radar.adc_full_scale = 2;
default.radar.rx_paths = [1:22,24,23];
default.radar.wfs(1).chan_equal_dB = [-0.7 -0.7 -0.7 -0.2 -0.1 -2.3 0.2 -0.1 -4.2 -4.5 -3.1 -4.6 -1.1 -1.5 -0.9 -1.7 0 -0.6 0.2 3.4 0.7 -2.1 0.2 0.9];
default.radar.wfs(1).chan_equal_deg = [-168.6 -114.1 -5.7 9 30 24.1 -144.3 -137.7 113.1 64.8 124.3 133.7 108 138.1 71.6 102.9 -95.8 -127.2 -143.1 -139.6 -128.4 -172.3 -158.5 -152.4];
default.radar.wfs(1).chan_equal_Tsys = [0.3 0.7 0 0.2 0.1 0.2 0.2 0.3 -31.8 -32.1 -31.7 -31.6 -31.2 -30.8 -31.5 -31.1 -4.3 -4.5 -4.6 -4.5 -4.5 -4.6 -4.5 -4.5]/1e9;
default.radar.noise_figure = 2;
default.radar.rx_gain = 10^(48/20);
default.radar.adc_SNR_dB = 59;
