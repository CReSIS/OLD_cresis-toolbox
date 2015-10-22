% script fix_records_for_leap_second
%
% The July 2012 leap second was not accounted for. This requires the
% following fixes to be done:
%
% 1. All GPS file should be remade with make_gps (do this step first)
%    THIS FUNCTION DOES NOT DO THIS
%
% 2a. For GPS files that used UTC time, the one second error will affect
%    the radar data AND the gps data in the same way so the only correction
%    to the radar data is to add one second to the GPS time.
%    THIS FUNCTION DOES THIS
%      accum: not implemented yet
%      fmcw: completed for records, data not implemented yet
%      rds: not implemented yet
%     
% 2b. For GPS files that used GPS time, the data should be completely
%    reprocessed since motion compensation may be off by one second.
%    The UTC time offset in the vectors sheet is likely to have a one
%    second error since if time synchronization was done correctly there
%    would have been a need to add one second to the UTC time to make
%    the radar data line up with the GPS data.  Therefore the UTC time
%    offset should be reverified.
%    THIS FUNCTION DOES NOT DO THIS
%
% Seasons affected and need to be fixed by 2a:
%   2012 Antarctica DC8, 2013 Greenland P3, Conway Accum???,
%   2013 Antarctica P3, 2014 Greenland P3
% Seasons affected and need to be fixed by 2b:
%   2013 Antarctica Basler*, 2013 Antarctica Ground*
%
% Author: John Paden

%% Settings

% params = read_param_xls(ct_filename_param('kuband_param_2012_Antarctica_DC8.xls'));
% params = read_param_xls(ct_filename_param('rds_param_2012_Antarctica_DC8.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2012_Antarctica_DC8.xls'));
% params = read_param_xls(ct_filename_param('accum_param_2013_Greenland_P3.xls'));
% params = read_param_xls(ct_filename_param('kuband_param_2013_Greenland_P3.xls'));
% params = read_param_xls(ct_filename_param('rds_param_2013_Greenland_P3.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2013_Greenland_P3.xls'));
params = read_param_xls(ct_filename_param('kuband_param_2013_Antarctica_Basler.xls'));
% params = read_param_xls(ct_filename_param('rds_param_2013_Antarctica_Basler.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2013_Antarctica_Basler.xls'));

update_records_en = false;
update_data_files_en = true;
% data_types = {'CSARP_post/layerData'};
% data_types = {'CSARP_post/qlook'};
data_types = {'qlook'};

%% Automated Section
if update_records_en
  for param_idx = 1:length(params)
    param = params(param_idx);
    if ~isfield(param.cmd,'generic') || ischar(param.cmd.generic) || ~param.cmd.generic
      continue;
    end
    fprintf('Updating records leap second %s\n', param.day_seg)
    records_fn = ct_filename_support(param,'','records');
    fprintf('  %s\n', records_fn);
    records = load(records_fn);
    [output_dir,radar_type] = ct_output_dir(param.radar_name);
    gps_fn = ct_filename_support(param,'','gps',1);
    fprintf('  %s\n', gps_fn);
    gps = load(gps_fn);
    if ~isfield(gps,'sw_version') || gps.sw_version.rev < 2332
      warning('GPS file has old software version, probably has wrong UTC leap second. Run make_gps before continuing.');
      keyboard
    end
    if records.param_records.sw_version.rev >= 2328
      warning('This records file was created with a sw version that had the UTC leap second correction. Are you sure you want to continue?');
      keyboard
    end
    
    if strcmpi(radar_type,'fmcw')
      utc_time_sod = double(records.raw.seconds) + double(records.raw.fraction)/param.radar.fs;
      
      %% Check for seconds of day roll over and unwrap (assume jump backward
      % of more than 23 hours is a roll over)
      wrap_idxs = find(abs(diff(utc_time_sod) + 86400) < 3600);
      for wrap_idx = wrap_idxs
        utc_time_sod(wrap_idx+1:end) = utc_time_sod(wrap_idx+1:end) + 86400;
      end
      
      %% Apply GPS sync correction to radar time
      utc_time_sod = utc_time_sod + param.vectors.gps.time_offset;
      
      %% Determine absolute radar time and convert from UTC to GPS
      year = str2double(param.day_seg(1:4));
      month = str2double(param.day_seg(5:6));
      day = str2double(param.day_seg(7:8));
      radar_gps_time = datenum_to_epoch(datenum(year,month,day,0,0,utc_time_sod)) + utc_leap_seconds(gps.gps_time(1));
    elseif strcmpi(output_dir,'rds')
      if strcmpi(param.radar_name,'mcords4')
        utc_time_sod = double(records.raw.seconds) + double(records.raw.fractions)/(param.radar.fs/4);
      else
        utc_time_sod = double(records.raw.seconds) + double(records.raw.fractions)/param.radar.fs;
      end
      
      %% Check for seconds of day roll over and unwrap (assume jump backward
      % of more than 23 hours is a roll over)
      wrap_idxs = find(abs(diff(utc_time_sod) + 86400) < 3600);
      for wrap_idx = wrap_idxs
        utc_time_sod(wrap_idx+1:end) = utc_time_sod(wrap_idx+1:end) + 86400;
      end
      
      %% Apply GPS sync correction to radar time
      utc_time_sod = utc_time_sod + param.vectors.gps.time_offset;
      
      %% Determine absolute radar time and convert from UTC to GPS
      year = str2double(param.day_seg(1:4));
      month = str2double(param.day_seg(5:6));
      day = str2double(param.day_seg(7:8));
      radar_gps_time = datenum_to_epoch(datenum(year,month,day,0,0,utc_time_sod)) + utc_leap_seconds(gps.gps_time(1));
    end
    
    find_utc_correction_note = regexpi(records.notes,'UTC Leap');
    if ~isempty(find_utc_correction_note)
      warning('UTC leap second correction is in the records.notes already. Are you sure you want to run again?');
      keyboard
    end
    
    records.notes = cat(2,records.notes,sprintf('\nUTC Leap Seconds corrected using %s\n', mfilename()));
    
    if strcmpi(param.season_name,'2013_Antarctica_Basler') ...
        && strcmpi(param.radar_name,'mcords4')
      time_correction = radar_gps_time(1) - records.gps_time(1);
      if abs(time_correction) > 0.2
        warning('Expecting there to be about a +0 second correction to apply to GPS time, but requires %f sec.', time_correction);
        keyboard
      end
    elseif ~strcmpi(output_dir,'accum')
      time_correction = radar_gps_time(1) - records.gps_time(1);
      if abs(time_correction - 1) > 0.2
        warning('Expecting there to be about a +1 second correction to apply to GPS time, but requires %f sec.', time_correction);
        %keyboard
      end
    else
      time_correction = mean(records.gps_time - records.raw.comp_time) - mean(gps.sync_gps_time-gps.comp_time)
      if abs(time_correction - 1) > 0.2
        warning('Expecting there to be about a +1 second correction to apply to GPS time, but requires %f sec.', time_correction);
        keyboard
      end
    end
    
    if strcmpi(param.season_name,'2013_Antarctica_Basler') ...
        && strcmpi(param.radar_name,'mcords4')
      records.param_records.vectors.gps.time_offset = records.param_records.vectors.gps.time_offset - 1;
      
      save(records_fn, '-append', '-struct', 'records', 'notes', 'param_records');
    else
      records.gps_time = records.gps_time + 1;
      
      save(records_fn, '-append', '-struct', 'records', 'notes', 'gps_time');
    end
  end
end

if update_data_files_en
  for param_idx = 1:length(params)
    param = params(param_idx);
    if ~isfield(param.cmd,'generic') || ischar(param.cmd.generic) || ~param.cmd.generic
      continue;
    end
    for data_types_idx = 1:length(data_types)
      data_type = data_types{data_types_idx};
      fprintf('Updating %s leap second %s\n', data_type, param.day_seg)
      data_fn_dir = ct_filename_out(param,data_type,'');
      fns = get_filenames(data_fn_dir,'Data_','','.mat');
      for fn_idx = 1:length(fns)
        fn = fns{fn_idx};
        fprintf('  %s (%s)\n', fn, datestr(now,'HH:MM:SS'));
        
        mdata = load(fn,'GPS_time','param_records');
        
        if ~isfield(mdata,'param_records')
          % This is a layer file
          
          if isfield(mdata,'utc_leap_second_correction')
            if 1
              warning('UTC leap second correction is in the data file already. Are you sure you want to run again?');
              keyboard
            else
              warning('UTC leap second correction is in the data file already. Skipping...');
              continue;
            end
          end
          mdata.utc_leap_second_correction = sprintf(', UTC Leap Seconds corrected using %s', mfilename());
          
          mdata.GPS_time = mdata.GPS_time + 1;
          
          save(fn, '-v6', '-append', '-struct', 'mdata', 'GPS_time', 'utc_leap_second_correction');
          
        else
          % This is a regular data file
          
          if ~isfield(mdata.param_records,'cmd')
            mdata.param_records.cmd.notes = '';
          end
          
          find_utc_correction_note = regexpi(mdata.param_records.cmd.notes,'UTC Leap');
          if ~isempty(find_utc_correction_note)
            if 0
              warning('UTC leap second correction is in the data file already. Are you sure you want to run again?');
              keyboard
            else
              warning('UTC leap second correction is in the data file already. Skipping...');
              continue;
            end
          end
          
          mdata.param_records.cmd.notes = cat(2,mdata.param_records.cmd.notes, ...
            sprintf(', UTC Leap Seconds corrected using %s', mfilename()));
          
          mdata.GPS_time = mdata.GPS_time + 1;
          
          save(fn, '-v6', '-append', '-struct', 'mdata', 'GPS_time', 'param_records');
        end
      end
    end
  end
end

return;
