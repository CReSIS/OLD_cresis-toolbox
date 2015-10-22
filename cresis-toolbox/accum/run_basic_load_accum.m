% script run_basic_load_accum
%
% Demonstrates how to use basic_load_accum and then process the data.
%
% Author: John Paden
%
% See also: basic_load_accum

% =======================================================================
% User Settings
% =======================================================================

param_file_version = '1.0';
param.radar_name = 'accum';
param.season_name = '2011_Greenland_P3';
param.day_seg = '20110412_01';

base_path = '/mnt/paden/20110418-accum/';

fn = get_filename(base_path,'accum','','0050.dat');

% Processing params(idx)eters
param.qlook.tukey = 0.5;
param.qlook.band_window_func = @hanning;
param.qlook.td_window_func = [];
param.qlook.window_func = @blackman;
param.qlook.sw_presums = 10;
param.qlook.incoh_ave = 10;
param.qlook.decimate = 5;

% Output parameters
param.qlook.out.en = false;
param.qlook.out.dir = '';
param.qlook.plot.en = true;
param.qlook.elev_comp.en = false;

param.qlook.gps.fn = '';
param.qlook.gps.time_offset = -1;

% Radar configuration
param.radar.fs = 1e9/8;
param.radar.f0 = 900e6;
param.radar.f_step = -20e6;
param.radar.BW = 50e6;
param.radar.Tpd = 2.048e-6;
param.radar.fLO = param.radar.f0 + 6.25e6;
param.radar.chan_equal = 10.^([0	0.884695377	2.969666953	4.609239012	4.270122918	3.519328706	4.534573365	4.361630813	3.883828894	5.188441292	5.675117442	4.882502887	4.187331947	4.98803758	5.700863022	4.737899962]/20);
param.radar.adc_bits = 14;
param.radar.Vpp_scale = 2;
param.radar.hw_presums = 8;
param.radar.bit_shifts = 5;
param.radar.wfs = 1:16;

param_override = [];

% =======================================================================
% Load data
% =======================================================================

[data,time] = qlook_accum_task(fn,param);

imagesc(lp(data));
colormap(1-gray(256));

return;
