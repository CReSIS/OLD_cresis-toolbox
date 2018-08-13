% convert_icards_to_new_format
%
% Function for converting old posted files to the new format (Aug 25, 2011).
%
% Special Notes for 1996_Greenland_P3:
%   load('Data_19960520_01_010.mat');
%   Latitude(500) = mean(Latitude([499 501]));
%   Longitude([490:491 500]) = interp1([1:489 492:499 501:length(Longitude)],Longitude([1:489 492:499 501:length(Longitude)]),[490:491 500]);
%   save('-append','Data_19960520_01_010.mat','Latitude','Longitude');
%
% Special Notes for 1999_Greenland_P3:
%   convert_icards_to_new_format_19990514
%   AND
%   rm /cresis/scratch2/mdce/icards/1999_Greenland_P3/CSARP_layerData/19990507_01/Data_19990507_01_008.mat
%   rm /cresis/scratch2/mdce/icards/1999_Greenland_P3/CSARP_standard/19990507_01/Data_19990507_01_008.mat
%
% Special Notes for 2002_Antarctica_P3chile:
%   frms = {};
%   frms{end+1} = '20021126_01_004';
%   frms{end+1} = '20021128_01_002';
%   frms{end+1} = '20021128_01_004';
%   frms{end+1} = '20021128_01_008';
%   frms{end+1} = '20021204_01_002';
%   frms{end+1} = '20021204_01_003';
%   frms{end+1} = '20021204_01_007';
%   frms{end+1} = '20021204_01_010';
%   frms{end+1} = '20021206_01_001';
%   frms{end+1} = '20021210_01_010';
%   frms{end+1} = '20021210_01_012';
%   frms{end+1} = '20021212_01_005';
%   frms{end+1} = '20021212_01_007';
%   base_path = '/cresis/scratch2/mdce/icards/2002_Antarctica_P3chile/CSARP_standard/';
%   for frm_idx = 1:length(frms)
%     frm = frms{frm_idx};
%     load(fullfile(base_path,frm(1:11),['Data_' frm '.mat']));
% bad_mask = abs(medfilt1(Latitude) - Latitude) > 0.1;
% figure(1); clf;
% plot(Latitude);
% hold on;
% Latitude(bad_mask) = interp1(find(~bad_mask),Latitude(find(~bad_mask)),find(bad_mask),'linear','extrap');
% plot(Latitude);
% hold off;
%     keyboard;
%   end
%
% Function is currently modified to ingest 2002chile
%
% Four ground based missions for ICARDS: 1999ngrip, 2003ngrip,
% 2004south_pole, 2006flade_raw
%   These are not done yet.
% 1999ngrip: not posted
% 2003ngrip: not posted
% 2004south_pole: not posted
% 2006flade_raw:
%   /cresis/data1/ICARDS/2006flade_raw/MAY25_06/
%   MAY25_06.0
%   MAY25_1.0
%   /cresis/data1/ICARDS/2006flade_raw/MAY26_06/
%   MAY26_1.0 
%   MAY26_2.0
%   /cresis/data1/ICARDS/2006flade_raw/MAY26__0/
%   MAY26_06.0 
%   /cresis/data1/ICARDS/2006flade_raw/MAY30_06/
%   MAY30_06.0
%   /cresis/data1/ICARDS/2006flade_raw/MAY31_06/
%   MAY31_06.0
%   /cresis/web/cresis_data/datafiles/greenland/2006_land-based/
%   Flight_Lines  MAT  May25_06  May26_06  May30_06  May31_06  PDF  TXT
%
% Author: John Paden, Shashanka Jagarlapudi

param.radar_name = 'icards';
param.year = 2002;
% param.season_name = sprintf('%04d_Greenland_P3',param.year);
param.season_name = sprintf('%04d_Antarctica_P3chile',param.year);

% data_base = sprintf('/cresis/web/cresis_data/datafiles/greenland/%04d/mat',param.year);
data_base = sprintf('/cresis/web/cresis_data/datafiles/Antarctica/%04d/MAT',param.year);

% thick_base = sprintf('/cresis/web/cresis_data/datafiles/greenland/%04d/thick',param.year);
thick_base = sprintf('/cresis/web/cresis_data/datafiles/Antarctica/%04d/TXT',param.year);

% raw_base = sprintf('/cresis/data1/ICARDS/%04d/',param.year);
raw_base = sprintf('/cresis/data1/ICARDS/%04dchile/data/',param.year);

param.out_dir = 'standard';

medfilt1_en = 1;

gps_med_filt_en = true;
if param.year == 1993 || param.year == 1995
  gps_checksum_en = false; % Only 1993,1995 must have this be false
else
  gps_checksum_en = true; % Only 1993,1995 must have this be false
end

% Old ICARDS .mat files have wrong sign for longitude usually, so setting
% this to -1 corrects it
param.longitude_sign = -1;

file_size_threshold = 4e3;
fs = 18.75e6;

% =======================================================================
% =======================================================================
% Automated Section
% =======================================================================
% =======================================================================

convert_icards_to_new_format_tstart = tic;

physical_constants;

data_dirs = get_filenames(data_base,'','','',struct('type','d'));
base_idx = strmatch(data_base,data_dirs,'exact');
data_dirs = data_dirs([1:base_idx-1, base_idx+1:end]);

for data_dir_idx = 1:length(data_dirs)
  data_dir = data_dirs{data_dir_idx};
  
  % Convert each directory
  fprintf('Converting directory %s\n', data_dir);
  
  [data_dir_dir data_dir_name] = fileparts(data_dir);
  
  try
%     data_date = datenum(sprintf('%s-%s-%04d',data_dir_name(4:5),data_dir_name(1:3),param.year));
    data_date = datenum(data_dir_name,'yyyymmdd');
  catch
    fprintf('  datenum conversion failed\n');
    continue
  end
  [param.year param.month param.day] = datevec(data_date);
  
  seg_num = 1;
  param.day_seg = sprintf('%04d%02d%02d_%02d',param.year,param.month,param.day,seg_num);
  
  out_dir = ct_filename_out(param, ...
    param.out_dir, 'CSARP_standard');
  
  % =======================================================================
  % Load txt file (to get GPS time)
  % =======================================================================
%   txt_fn = fullfile(thick_base,sprintf('%s_%02d.txt',data_dir_name,mod(param.year,100)));
  txt_fn = fullfile(thick_base,sprintf('%s.txt',data_dir_name));
  fid = fopen(txt_fn,'r');
  clear txt;
  if fid > 1
    txt.data = textscan(fid,'%f %f %f %f %f %s','Headerlines',1);
    fclose(fid);
    txt.lat = round(txt.data{1}.'*1e6)/1e6;
    txt.lon = round(txt.data{2}.'*1e6)/1e6;
    txt.utc_time_sod = round(txt.data{3}.'*1e6)/1e6;
    txt.thick = round(txt.data{4}.'*1e6)/1e6;
    txt.height = round(txt.data{5}.'*1e6)/1e6;
  end
  
  % =======================================================================
  % Load raw data files (to get raw GPS)
  %   Raw data files have the NMEA string stored in them
  %   The NMEA string can have a lot of errors, we just throw out
  %   strings with bad errors.
  % =======================================================================
  clear gps;
  gps.gps_time = [];
  gps.lat = [];
  gps.lon = [];
  gps.elev = [];
  raw_dir_day = datestr(data_date,'dd');
  if raw_dir_day(1) == '0'
    raw_dir_day = raw_dir_day(2);
  end
  raw_dir = fullfile(raw_base,[lower(datestr(data_date,'mmm')) raw_dir_day]);
  fprintf('  Loading GPS data from raw data: %s\n', raw_dir);
  fns = get_filenames(raw_dir,lower(datestr(data_date,'mmm')),raw_dir_day,'.[0-9][0-9][0-9]',struct('exact',1));
  gps_rline = 0;
  for fn_idx = 1:length(fns)
    fn = fns{fn_idx};
    [fn_dir,fn_name,fn_ext] = fileparts(fn);
    fn_name = [fn_name fn_ext];
    file_info = dir(fn);
    fprintf('  %s (%d), (%.1f sec)\n', fn_name, file_info.bytes, toc(convert_icards_to_new_format_tstart));
    gps_data = char(icards_get_data(fn,4).');
    gps_offset_idx = length(gps.lat);
    gps.gps_time = cat(2,gps.gps_time,NaN*zeros(1,size(gps_data,1)));
    gps.lat = cat(2,gps.lat,NaN*zeros(1,size(gps_data,1)));
    gps.lon = cat(2,gps.lon,NaN*zeros(1,size(gps_data,1)));
    gps.elev = cat(2,gps.elev,NaN*zeros(1,size(gps_data,1)));
    for rline = 1:size(gps_data,1)
      % Parse each line of GPS input, gps_data(rline,:)
      C = textscan(gps_data(rline,:),'%s','Delimiter',',');
      gps_rline = gps_rline+1;
      if gps_checksum_en
        try
          checksum = hex2dec(gps_data(rline,end-1:end));
        catch
          continue;
        end
        if bin2dec(char(mod(sum(dec2bin(gps_data(rline,2:end-3))-'0'),2)+'0')) ~= checksum
          continue;
        end
      end
      if length(C) < 1 || length(C{1}) < 10
        continue;
      end
      if ~strcmp(C{1}{1},'$GPGGA')
        continue;
      end
      if length(C{1}{2}) ~= 9 || C{1}{2}(1) >= '3' || C{1}{2}(7) ~= '.' || length(find(C{1}{2}=='.')) > 1
        continue;
      end
      if length(C{1}{3}) < 9 || length(C{1}{3}) > 10
        continue;
      end
      if length(C{1}{4}) ~= 1
        continue;
      end
      if length(C{1}{5}) < 10 || length(C{1}{5}) > 11
        continue;
      end
      if length(C{1}{6}) ~= 1
        continue;
      end
      if length(C{1}{10}) ~= 6 && length(C{1}{10}) ~= 9
        continue;
      end
      hour = str2double(C{1}{2}(1:2));
      minute = str2double(C{1}{2}(3:4));
      second = str2double(C{1}{2}(5:end));
      if hour > 23 || minute > 59 || second > 59.99
        continue;
      end
      gps_time_epoch = datenum_to_epoch(datenum(param.year,param.month, ...
        param.day,hour,minute,second));
      gps.gps_time(gps_rline) = gps_time_epoch;
      gps.lat(gps_rline) = str2double(C{1}{3}(1:2)) + str2double(C{1}{3}(3:end))/60;
      if C{1}{4} ~= 'N'
        gps.lat(gps_rline) = -gps.lat(gps_rline);
      end
      gps.lon(gps_rline) = str2double(C{1}{5}(1:3)) + str2double(C{1}{5}(4:end))/60;
      if C{1}{6} ~= 'E'
        gps.lon(gps_rline) = -gps.lon(gps_rline);
      end
      gps.elev(gps_rline) = str2double(C{1}{10});
    end
  end
  
  save('/cresis/scratch1/paden/mdce_tmp/convert_icards');
  good_mask = ~isnan(gps.gps_time);
  gps.gps_time = gps.gps_time(good_mask);
  gps.lat = gps.lat(good_mask);
  gps.lon = gps.lon(good_mask);
  gps.elev = gps.elev(good_mask);
  
  if gps_med_filt_en
    gps.gps_time = medfilt1(gps.gps_time,5);
    gps.lat = medfilt1(gps.lat,5);
    gps.lon = medfilt1(gps.lon,5);
    gps.elev = medfilt1(gps.elev,5);
  end
  
  % Force GPS time to be monotonically increasing
  good_mask = ones(size(gps.gps_time));
  cur_time = gps.gps_time(1);
  for gps_rline = 2:length(gps.gps_time)
    if gps.gps_time(gps_rline) <= cur_time
      good_mask(gps_rline) = 0;
    else
      cur_time = gps.gps_time(gps_rline);
    end
  end
  good_mask = logical(good_mask);
  gps.gps_time = gps.gps_time(good_mask);
  gps.lat = gps.lat(good_mask);
  gps.lon = gps.lon(good_mask);
  gps.elev = gps.elev(good_mask);
  
  % NMEA strings are supposed to be UTC time, so convert to GPS time
  gps.gps_time = gps.gps_time + utc_leap_seconds(gps.gps_time(1));
  
  % Find ECEF fixed for each position with elevation set to zero
  [gps.x,gps.y,gps.z] = geodetic2ecef(gps.lat/180*pi,gps.lon/180*pi, ...
    0*gps.elev,WGS84.ellipsoid);
  
  figure(1); clf; plot(gps.gps_time);
  figure(2); clf; plot(gps.lat);
  drawnow;
  
  % =======================================================================
  % Load output data (.mat) files
  % =======================================================================
  fns = get_filenames(data_dir,'','','.mat');
  
  frm = 0;
  for fn_idx = 1:length(fns)
    fn = fns{fn_idx};
    [fn_dir,fn_name] = fileparts(fn);
    file_info = dir(fn);
    fprintf('  %s (%d)\n', fn_name, file_info.bytes);
    
    % Remove files that are too small
    if file_info.bytes < file_size_threshold
      fprintf('    Too small\n');
      continue;
    end
    
    % Remove files with non-conforming name
%     if length(fn_name) ~= 25
    if length(fn_name) ~= 15 + 2*length(raw_dir_day)
      fprintf('    Wrong filename length\n');
      continue;
    end
    
    frm = frm + 1;
    
    % Load file
    tmp = load(fn);
    
    % Correct longitude sign if necessary
    tmp.longitude = param.longitude_sign * tmp.longitude;
    
    Data = tmp.A;
    dt = 1/fs;
    Nt = size(Data,1);
    Time = dt*(0:Nt-1).';
    Depth = Time * c/2;
    Surface = interp1(0:Nt-1,Time,tmp.top);
    Bottom = interp1(0:Nt-1,Time,tmp.bot);
    
    if fn_idx == 23 && param.year == 1999 && param.month == 5 && param.day == 26
      txt.utc_time_sod = medfilt1(txt.utc_time_sod,7);
    end
    
    if exist('txt','var')
      % Output thickness (.txt) file exists, since these files are
      % (supposed to be) exactly the same at the .mat files we can
      % use the .txt file to get GPS_time
      comp_lat = round(tmp.latitude*1e4)/1e4;
      comp_lon = round(tmp.longitude*1e4)/1e4;
      
      mean_error = zeros(length(txt.lat)-length(tmp.latitude)+1,1);
      for offset = 1:length(txt.lat)-length(tmp.latitude)+1
        mean_error(offset) = sum(abs(txt.lat(offset + (0:length(comp_lat)-1)) ...
          - comp_lat));
      end
      [min_mean_error offset] = min(mean_error);
      utc_time_sod = txt.utc_time_sod(offset + (0:length(comp_lat)-1));
      utc_time = datenum_to_epoch(datenum(param.year,param.month, ...
        param.day,0,0,utc_time_sod));
      GPS_time = utc_time + utc_leap_seconds(utc_time(1));
    else
      % Output thickness (.txt) file does not exist, use raw data to
      % find GPS_time
      fprintf('    No thickness file! Using less accurate raw file sync\n');
      [tmp.x,tmp.y,tmp.z] = geodetic2ecef(tmp.latitude/180*pi,tmp.longitude/180*pi, ...
        zeros(size(tmp.latitude)),WGS84.ellipsoid);
      clear min_idx min_dist;
      for rline = 1:size(Data,2)
        % Calculate the distance to every gps point recorded, set the GPS
        % time to the one that is closest
        dist = sqrt((tmp.x(rline)-gps.x).^2 + (tmp.y(rline)-gps.y).^2 ...
          + (tmp.z(rline)-gps.z).^2);
        [min_dist(rline) min_idx(rline)] = min(dist);
      end
      GPS_time = gps.gps_time(medfilt1(min_idx,7));
    end
    
    Latitude = tmp.latitude;
    Longitude = tmp.longitude;
    if medfilt1_en
      bad_mask = abs(medfilt1(GPS_time,3) - GPS_time) > 10*median(diff(GPS_time));
      if any(bad_mask)
        find(bad_mask)
        figure(1); clf;
        plot(GPS_time);
        GPS_time(bad_mask) = interp1(find(~bad_mask),GPS_time(find(~bad_mask)),find(bad_mask),'linear','extrap');
        hold on;
        plot(GPS_time,'r');
        hold off;
        keyboard
      end
    end
    Elevation = interp1(gps.gps_time,gps.elev,GPS_time,'linear','extrap');
    
    % Force GPS times to be monotonically increasing
    for rline = 2:length(GPS_time)
      if GPS_time(rline) <= GPS_time(rline-1)
        GPS_time(rline) = GPS_time(rline-1) + 1e-4;
      end
    end
    
    out_fn_name = sprintf('Data_%s_%03d.mat',param.day_seg,frm);
    out_fn = fullfile(out_dir,out_fn_name);
    if ~exist(out_dir,'dir')
      mkdir(out_dir);
    end
    fprintf('    Saving %s\n', out_fn);
    save(out_fn,'Data','Time','Depth','Latitude','Longitude','Elevation','GPS_time','Surface','Bottom');
  end
  keyboard
  
end

return;
