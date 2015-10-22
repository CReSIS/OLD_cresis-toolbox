% update_sar_files
%
% Used for updating sar files (fk_data) with new GPS information,
% GPS time offset, or layer information.
%
% Author: John Paden

%fn_dir = 'V:/output/mcords2/2011_Greenland_P3/CSARP_out/20110507_02';
fn_dir = '/cresis/scratch2/mdce/mcords2/2011_Greenland_P3/CSARP_out/20110413_03';
new_time_offset = -1;

% gps_fn = 'V:/csarp_support/gps/2011_Greenland_P3/gps_20110507.mat';
gps_fn = '/cresis/scratch1/mdce/csarp_support/gps/2011_Greenland_P3/gps_20110413.mat';

param.frms = [14:16 18:20]; % Leave empty to do all

param.radar_name = 'mcords2';
param.season_name = '2011_Greenland_P3';

param.gps.update_en = false;

param.layers.update_en = true;
param.layers.path = '';

param.save_en = true;

% ======================================================================
% ======================================================================

update_sar_files_tstart = tic;

fns = get_filenames(fn_dir,'fk_data','','',struct('type','d'));

gps = load(gps_fn);

if param.layers.update_en
  physical_constants;
  
  % Load all the layer data
  layer_path = ct_filename_out(param,param.layers.path,'CSARP_layerData',1);
  layer_fns = get_filenames(layer_path,'Data','','.mat',struct('recursive',1));
  
  fprintf('Loading layer files (%.1f sec)\n', toc(update_sar_files_tstart));
  surface = [];
  bottom = [];
  elev = [];
  GPS_time = [];
  for file_idx = 1:length(layer_fns)
    tmp = load(layer_fns{file_idx});
    surface = [surface tmp.layerData{1}.value{2}.data];
    bottom = [bottom tmp.layerData{2}.value{2}.data];
    elev = [elev tmp.Elevation];
    GPS_time = [GPS_time tmp.GPS_time];
    if length(tmp.GPS_time) ~= length(tmp.layerData{1}.value{2}.data)
      fprintf('GPS_time dimensions do not match surface dimensions (ERROR!)\n')
      fprintf('  %s\n', layer_fns{file_idx});
      return;
    end
  end
  
  % Load of layer files may not be in order of the time that the data
  % was collected, so we fix that here since we will be interpolating
  % using the time axis later.
  [GPS_time,sorting_idxs] = sort(GPS_time);
  surface = surface(sorting_idxs);
  bottom = bottom(sorting_idxs);
  elev = elev(sorting_idxs);
end

if param.gps.update_en
  % Get the GPS seconds of day to sync to radar
  utc_time_datenum = epoch_to_datenum(gps.gps_time - utc_leap_seconds(gps.gps_time(1)));
  [year month day hour minute sec] = datevec(utc_time_datenum);
  UTC_sod = (day-day(1))*86400+hour*3600+minute*60+sec;  % GPS seconds of day
end

%
for file_idx = 1:length(fns)
  fn = fns{file_idx};
  
  % =======================================================================
  % Select only the frames files that are in the frames list
  % =======================================================================
  [fn_dir fn_name] = fileparts(fn);
  frm = str2double(fn_name(9:11));
  if isfield(param,'frms') && ~isempty(param.frms)
    if ~any(frm == param.frms)
      continue;
    end
  end
  fprintf('Updating frame %s, %d of %d (%.1f sec)\n', fn_name, ...
    file_idx, length(fns), toc(update_sar_files_tstart));
  
  fk_fns = get_filenames(fn,'','','*.mat');
  
  for fk_idx = 1:length(fk_fns)
    fk_fn = fk_fns{fk_idx};
    
    [fk_fn_dir fk_fn_name] = fileparts(fk_fn);
    fprintf('  Updating file %s, %d of %d (%.1f sec)\n', fk_fn_name, ...
      fk_idx, length(fk_fns), toc(update_sar_files_tstart));
    
    load(fk_fn);
    
    if param.gps.update_en
      delta_offset = new_time_offset - param_records.gps.time_offset;
      param_records.gps.time_offset = new_time_offset;
      records.gps_time = records.gps_time + delta_offset;
      records.lat = interp1(gps.gps_time, gps.lat, records.gps_time);
      records.lon = interp1(gps.gps_time, gps.lon, records.gps_time);
      records.elev = interp1(gps.gps_time, gps.elev, records.gps_time);
      if param.save_en
        save(fk_fn,'records','param_records','-APPEND');
      else
        fprintf('  Not saving information (TEST MODE)\n');
      end
    end
    
    if param.layers.update_en
      % Update surface and bottom variables from layer files
      % 1. Interpolate layer file variables onto record timestamps
      new_elev = interp1(GPS_time,elev,records.gps_time);
      new_surface = interp1(GPS_time,surface,records.gps_time);
      new_bottom = interp1(GPS_time,bottom,records.gps_time);
      % 2. Adjust for differences in measurement elevation due to
      %    motion compensation
      new_surface = (records.elev-(new_elev-new_surface*c/2))*2/c;
      new_bottom = (records.elev-(new_elev-new_bottom*c/2))*2/c;
      tmp = records;
      % 3. Update surface
      %         if isfield(records,'surface')
      %           records.surface(isfinite(new_surface)) ...
      %             = new_surface(isfinite(new_surface));
      %         else
      %           records.surface = surface;
      %         end
      % 3. Update bottom
      if isfield(records,'bottom')
        records.bottom(isfinite(new_bottom)) ...
          = new_bottom(isfinite(new_bottom));
      else
        records.bottom = new_bottom;
      end
      if param.save_en
        save(fk_fn,'records','-APPEND');
      else
        fprintf('  Not saving information (TEST MODE)\n');
      end
    end
  end
  
end

return;






