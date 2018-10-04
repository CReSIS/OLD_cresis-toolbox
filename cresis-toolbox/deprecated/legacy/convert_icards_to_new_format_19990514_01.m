% ================================================================
% User Settings
% ================================================================

param.radar_name = 'icards';
param.year = 1999;
param.month = 5;
param.day = 14;
param.season_name = sprintf('%04d_Greenland_P3',param.year);

raw_base = sprintf('/cresis/data1/ICARDS/%04d/',param.year);
data_dir_name = 'may14';
param.day_seg = '19990514_01';

param.out_dir = 'standard';

gps_med_filt_en = true;
if param.year == 1993 || param.year == 1995
  gps_checksum_en = false; % Only 1993,1995 must have this be false
else
  gps_checksum_en = true; % Only 1993,1995 must have this be false
end

file_size_threshold = 4e3;
fs = 18.75e6;
Mx = 10; % slow-time decimation

frms  = [23 24 25 26];
start = [71 74 77 80];
stop  = [73 76 79 83];

% ================================================================
% Automated section
% ================================================================

physical_constants;
clear gps;
raw_dir = fullfile(raw_base,data_dir_name);

out_dir = ct_filename_out(param, ...
  param.out_dir, 'CSARP_standard');

convert_icards_to_new_format_tstart = tic;

for frm_idx = 1:length(frms)
  frm = frms(frm_idx);
  data = [];
  gps.gps_time = [];
  gps.lat = [];
  gps.lon = [];
  gps.elev = [];
  gps_rline = 0;
  for file_idx = start(frm_idx) : stop(frm_idx)
    fn = get_filename(raw_dir,'','',sprintf('may14_99.%03.0f',file_idx));
    
    fprintf('Loading file %s (%.1f sec)\n', fn, toc(convert_icards_to_new_format_tstart));
    
    % ===============================================================
    % Load Radar Data
    % ===============================================================
    fields = icards_get_available(fn);
    if any(fields(:,1) == 3)
      coherent_en = true;
    else
      coherent_en = false;
    end
    data = cat(2,data,icards_get_data(fn,2)+j*icards_get_data(fn,3));
    
    % ===============================================================
    % Load GPS Data
    % ===============================================================
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
      if length(C{1}{2}) ~= 9 || C{1}{2}(1) >= '2' || C{1}{2}(7) ~= '.' || length(find(C{1}{2}=='.')) > 1
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
  % Remove DC bias
  %for rline = 1:size(data,2)
  %  data(:,rline) = data(:,rline) - mean(data(650:end,rline));
  %end
  
  % Incoherent averaging
  data = fir_dec(abs(data).^2,Mx);
  
  % GPS Cleanup
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
  good_mask_orig = good_mask;
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
  good_mask_final = good_mask_orig;
  good_mask_final(good_mask_orig) = good_mask;
  
  % NMEA strings are supposed to be UTC time, so convert to GPS time
  gps.gps_time = gps.gps_time + utc_leap_seconds(gps.gps_time(1));
  
  % Find ECEF fixed for each position with elevation set to zero
  [gps.x,gps.y,gps.z] = geodetic2ecef(gps.lat/180*pi,gps.lon/180*pi, ...
    0*gps.elev,WGS84.ellipsoid);
  
  GPS_time = fir_dec(interp1(find(good_mask_final),gps.gps_time,1:length(good_mask_final),'linear','extrap'),Mx);
  Latitude = fir_dec(interp1(find(good_mask_final),gps.lat,1:length(good_mask_final),'linear','extrap'),Mx);
  Longitude = fir_dec(interp1(find(good_mask_final),gps.lon,1:length(good_mask_final),'linear','extrap'),Mx);
  Elevation = fir_dec(interp1(find(good_mask_final),gps.elev,1:length(good_mask_final),'linear','extrap'),Mx);
  
  if 1
    figure(1); clf;
    imagesc(lp(data));
    colormap(1-gray(256));
    ylim([1 700]);
    
    figure(2); clf;
    plot(gps.gps_time);
    
    figure(3); clf;
    plot(gps.lat,gps.lon,'bx');
  end
  
  Data = data;
  dt = 1/fs;
  Nt = size(Data,1);
  Time = dt*(0:Nt-1).';
  Depth = Time * c/2;
  Surface = NaN*zeros(size(GPS_time));
  Bottom = NaN*zeros(size(GPS_time));
  
  out_fn_name = sprintf('Data_%s_%03d.mat',param.day_seg,frm);
  out_fn = fullfile(out_dir,out_fn_name);
  if ~exist(out_dir,'dir')
    mkdir(out_dir);
  end
  fprintf('    Saving %s\n', out_fn);
  save(out_fn,'Data','Time','Depth','Latitude','Longitude','Elevation','GPS_time','Surface','Bottom');
end

return;
