% script convert_mcrds_old_layer_files
%
% Takes old MCRDS layer files and overwrites the bottom data in the current MCRDS
% radar files with the old information. The script is currently set up to do
% 2009_Greenland_TO.
%
% Author: Kyle Purdon

return; % For safety (remove when you want to run this)

fprintf('==============================================\n');
fprintf('Starting layer conversion\n');
fprintf('==============================================\n');

old_dirs = {};
fast_time_corr = [];

old_dirs{end+1} = '/cresis/web/cresis_data/datafiles/greenland/2009/Jakobshavn/20090331/MAT/';
fast_time_corr(end+1) = 0e-6;
old_dirs{end+1} = '/cresis/web/cresis_data/datafiles/greenland/2009/Jakobshavn/20090401/MAT/';
fast_time_corr(end+1) = 0e-6;
old_dirs{end+1} = '/cresis/web/cresis_data/datafiles/greenland/2009/Jakobshavn/20090402/MAT/';
fast_time_corr(end+1) = 0e-6;
old_dirs{end+1} = '/cresis/web/cresis_data/datafiles/greenland/2009/Jakobshavn/20090403/MAT/';
fast_time_corr(end+1) = 0e-6;
old_dirs{end+1} = '/cresis/web/cresis_data/datafiles/greenland/2009/Helheim/20090409/MAT/';
fast_time_corr(end+1) = 0e-6;
old_dirs{end+1} = '/cresis/web/cresis_data/datafiles/greenland/2009/Helheim/20090411/MAT/';
fast_time_corr(end+1) = 0e-6;
old_dirs{end+1} = '/cresis/web/cresis_data/datafiles/greenland/2009/Helheim/20090414a/MAT/';
fast_time_corr(end+1) = 0e-6;
old_dirs{end+1} = '/cresis/web/cresis_data/datafiles/greenland/2009/Kangerdlugssuaq/20090414b/MAT/';
fast_time_corr(end+1) = 0e-6;
old_dirs{end+1} = '/cresis/web/cresis_data/datafiles/greenland/2009/Kangerdlugssuaq/20090422a/MAT/';
fast_time_corr(end+1) = 0e-6;
old_dirs{end+1} = '/cresis/web/cresis_data/datafiles/greenland/2009/Kangerdlugssuaq/20090422b/MAT/';
fast_time_corr(end+1) = 0e-6;
old_dirs{end+1} = '/cresis/web/cresis_data/datafiles/greenland/2009/Kangerdlugssuaq/20090423a/MAT/';
fast_time_corr(end+1) = 0e-6;

new_dirs ={};

new_dirs{end +1} = '/cresis/scratch2/mdce/mcrds/2009_Greenland_TO/CSARP_layerData/20090331_01/';
new_dirs{end +1} = '/cresis/scratch2/mdce/mcrds/2009_Greenland_TO/CSARP_layerData/20090402_01';
new_dirs{end +1} = '/cresis/scratch2/mdce/mcrds/2009_Greenland_TO/CSARP_layerData/20090406_01/';
new_dirs{end +1} = '/cresis/scratch2/mdce/mcrds/2009_Greenland_TO/CSARP_layerData/20090422_01/';
new_dirs{end +1} = '/cresis/scratch2/mdce/mcrds/2009_Greenland_TO/CSARP_layerData/20090422_02/';
new_dirs{end +1} = '/cresis/scratch2/mdce/mcrds/2009_Greenland_TO/CSARP_layerData/20090423_01/';
new_dirs{end +1} = '/cresis/scratch2/mdce/mcrds/2009_Greenland_TO/CSARP_layerData/20090423_02/';
new_dirs{end +1} = '/cresis/scratch2/mdce/mcrds/2009_Greenland_TO/CSARP_layerData/20090428_01/';

if 0
old_surface = [];
old_bottom = [];
old_UTC_time = [];
for dir_idx = 1:length(old_dirs)
  old_dir = old_dirs{dir_idx};
  fns = get_filenames(old_dir,'Data_2009','','.mat');
  for fn_idx = 1:length(fns)
    fn = fns{fn_idx};
    fprintf('Loading %s\n', fn);
    old = load(fn,'Surface','Bottom','UTC_time','Depth','Time');
    fprintf('  First time %.1f us\n', old.Time(1)*1e6);
    old.Time =  old.Time(1:length(old.Depth)) + (fast_time_corr(dir_idx) - old.Time(1));
    old_surface = cat(2,old_surface,interp1(old.Depth,old.Time,old.Surface));
    old_bottom = cat(2,old_bottom,interp1(old.Depth,old.Time,old.Bottom));
    old_UTC_time = cat(2,old_UTC_time,old.UTC_time);
  end
end
end

param.radar_name = 'mcrds';
param.season_name = '2009_Greenland_TO';

for dir_idx = 1:length(new_dirs)
  fns = get_filenames(new_dirs{dir_idx},'Data','','.mat');
  new_surface =[];
  new_bottom =[];
  new_GPS_time = [];
  for fn_idx = 1:length(fns)
    fn = fns{fn_idx};
    fprintf('Loading %s\n', fn);
    new = load(fn,'layerData','GPS_time');
    old = new;
%     new.layerData{1}.value{2}.data = interp1(old_UTC_time,old_surface,new.GPS_time - utc_leap_seconds(old_UTC_time(1)));
    new.layerData{2}.value{2}.data = interp1(old_UTC_time,old_bottom,new.GPS_time - utc_leap_seconds(old_UTC_time(1)));
    if 0
      [fn_dir fn_name] = fileparts(fn);
      param.day_seg = fn_name([6:16]);
      data_fn = fullfile(ct_filename_out(param,'standard'),[fn_name '.mat']);
      tmp = load(data_fn);
      new_GPS_time = linspace(tmp.GPS_time(1),tmp.GPS_time(end),length(tmp.GPS_time));
      new_data = interp1(tmp.GPS_time,tmp.Data.',new_GPS_time).';
      new_surface = interp1(tmp.GPS_time, new.layerData{1}.value{2}.data, new_GPS_time);
      new_bottom = interp1(tmp.GPS_time, new.layerData{2}.value{2}.data, new_GPS_time);
      figure(1); clf;
      imagesc([],tmp.Time,lp(new_data));
      hold on;
      plot(new_surface-2.5e-6,'k.')
      plot(new_bottom-2.5e-6,'k.')
      hold off;
      ylim([-1e-6 30e-6]);
      keyboard
    end
    save(fn,'-struct','new','layerData');
  end
end









