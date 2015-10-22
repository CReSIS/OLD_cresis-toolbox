function [lat,lon,elev] = read_ramp_pass(ramp_fn)
% [lat,lon,elev] = read_ramp_pass(ramp_fn)
%
% Load ramp pass data file from ATM/John Sonntag
%
% Example:
%  ramp_fn = fullfile(gRadar.data_support_path,'ATM_ramp_passes','121005_truk_elmirage_anthtrem_nopark.txt');
%  [lat,lon,elev] = read_ramp_pass(ramp_fn);
%  scatter(lon,lat,[],elev,'.');
%
%  ramp_fn = fullfile(gRadar.data_support_path,'ATM_ramp_passes','121014_truk_L12_puntaarenas_anthtrem_nopark.txt');
%  [lat,lon,elev] = read_ramp_pass(ramp_fn);
%  scatter(lon,lat,[],elev,'.');
%
%  ramp_fn = fullfile(gRadar.data_support_path,'ATM_ramp_passes','130123_truk_l12_N159_anthtrem_nopark.txt');
%  [lat,lon,elev] = read_ramp_pass(ramp_fn);
%  scatter(lon,lat,[],elev,'.');
%
%  ramp_fn = fullfile(gRadar.data_support_path,'ATM_ramp_passes','130407_truk_l12_kangramp_anthtrem_nopark.txt');
%  [lat,lon,elev] = read_ramp_pass(ramp_fn);
%  scatter(lon,lat,[],elev,'.');
%
%  ramp_fn = fullfile(gRadar.data_support_path,'ATM_ramp_passes','130428_truk_l12_thuleramp_anthtrem_nopark.txt');
%  [lat,lon,elev] = read_ramp_pass(ramp_fn);
%  scatter(lon,lat,[],elev,'.');
%
%  ramp_fn = fullfile(gRadar.data_support_path,'ATM_ramp_passes','131114_atv_l12_jgs_itrf08_23dec13_blackislandtemp_anthtrem.txt');
%  [lat,lon,elev] = read_ramp_pass(ramp_fn);
%  scatter(lon,lat,[],elev,'.');
%
%  ramp_fn = fullfile(gRadar.data_support_path,'ATM_ramp_passes','131117_truk_l12_jgs_itrf08_30dec13_thel_anthtrem.txt');
%  [lat,lon,elev] = read_ramp_pass(ramp_fn);
%  scatter(lon,lat,[],elev,'.');
%
%  ramp_fn = fullfile(gRadar.data_support_path,'ATM_ramp_passes','131126_truk_l12_jgs_itrf08_30dec13_thel_anthtrem.txt');
%  [lat,lon,elev] = read_ramp_pass(ramp_fn);
%  scatter(lon,lat,[],elev,'.');
%
% Author: John Paden

fid = fopen(ramp_fn);
A = textscan(fid,'%f %f %f', 'Headerlines',1);
fclose(fid);
[lat,lon,elev] = deal(A{:});

end
