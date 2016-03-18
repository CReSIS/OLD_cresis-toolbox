% script deconv.remove_bad_wfs

%% User Settings
clear param;
param.radar_name = 'snow';
% param.season_name = '2009_Greenland_P3';
% param.day_seg = '20090323_09';
% param.season_name = '2010_Greenland_DC8';
% param.day_seg = '20100323_09';
% param.season_name = '2011_Greenland_P3';
% param.day_seg = '20110415_01';
% param.season_name = '2012_Greenland_P3';
% param.day_seg = '20120315_03';
param.season_name = 'SEASON';
param.day_seg = 'DAY_SEG';

bad_wfs = [4];

spec_file_input_type = 'noise';

%% Automated Section

fprintf('Loading %s\n', param.day_seg);
fn_dir = fileparts(ct_filename_out(param,spec_file_input_type, ''));
    
fn = fullfile(fn_dir,sprintf('deconv_%s.mat', param.day_seg));
    
final = load(fn);

final.deconv_gps_time(bad_wfs) = NaN;

fprintf('Type dbcont if you are sure you want to continue.\n');
keyboard
save(fn,'-struct','final');
