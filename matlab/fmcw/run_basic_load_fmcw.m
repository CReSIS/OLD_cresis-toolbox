% script run_basic_load_fmcw
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
param.radar_name = 'kuband';
param.season_name = '2011_Greenland_P3';
param.day_seg = '20110512_01';

base_path = '/mnt/paden/20110418-accum/';
base_path = 'd:\kuband\';

fn = get_filename(base_path,'kuband','','0290.dat');

param.qlook.gps.en = 0;
param.qlook.gps.fn = '';
param.qlook.gps.time_offset = 1;

param.qlook.out.en = 0;
param.qlook.out.dir = '';

param.qlook.window_func = @hanning;
param.qlook.st_window_func = [];
param.qlook.sw_presums = 4;
param.qlook.incoh_ave = [3 10];
param.qlook.decimate = 6;
param.qlook.plot.en = 2;
param.qlook.elev_comp.en = 0;
param.qlook.good_rbins = [751 15750];
param.qlook.nyquist_zone = 1;
filter_func = @() fir1(64,0.2,'high');
[param.qlook.B_st_filt,param.qlook.A_st_filt] = filter_func();
param.qlook.detrend_poly_order = 2;
param.qlook.surf.min_bin = 0;
param.qlook.surf.search_rng = [-120:120];

param.radar.prf = 2000;
param.radar.fs = 1e9/16;
param.radar.f0 = 232e6;
param.radar.f1 = 295e6;
param.radar.fmult = 56;
param.radar.fLO = 0e9;
param.radar.Tpd = 250e-6;
param.radar.adc_bits = 14;
param.radar.Vpp_scale = 2;
param.radar.td = 0;

param_override = [];

% =======================================================================
% Load data
% =======================================================================

physical_constants;

[data,time] = qlook_fmcw_task(fn,param);

range = c/2*time;
imagesc([], range, lp(data));
colormap(1-gray(256));

return;
