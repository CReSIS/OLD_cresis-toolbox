function [param,defaults] = default_radar_params_2016_Greenland_P3_mcords
% [param,defaults] = default_radar_params_2016_Greenland_P3_mcords
%
% Creates base "param" struct
% Creates defaults cell array
%
% Author: John Paden

param.season_name = '2016_Greenland_P3';
param.radar_name = 'mcords5';

default.radar.fs = 125e6;
default.radar.adc_bits = 14;
default.radar.adc_full_scale = 2;
default.radar.rx_paths = [1:15];
default.radar.wfs(1).chan_equal_dB = zeros(1,15);
default.radar.wfs(1).chan_equal_deg = zeros(1,15);
default.radar.wfs(1).chan_equal_Tsys = zeros(1,15);
default.radar.noise_figure = 2;
default.radar.rx_gain = 51.5;
default.radar.adc_SNR_dB = 70;
default.radar.Tadc_adjust = 0;

default.noise_50ohm = zeros(1,15);

default.Pt = (7*300) * sum(chebwin(7,30).^2)/7;

default.Gt = 7*4;
default.Ae = 2*0.468 * 0.468;

default.system_loss_dB = 10.^(-5.88/10);

defaults{1} = default;

return;
