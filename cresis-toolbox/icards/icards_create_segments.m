

% script create_segments_icards.m
%
% Creates the vectors worksheet for the param spreadsheet.
%
% Author: Kelly Rodriguez, John Paden

%to see extrapolated data using epri use may 14, 1991
%find date which non chronological data 
%% User Settings
param=[];
date=19970521;% 8 digits date                                %change needed 
date_string=num2str(date);
[month1,month2,day1,day2]=icards_monthANDday(date);
if ispc
    base_dir = strcat('Z:\ICARDS\',data_string(1:4),'\');
    adc_folder_name=strcat(month2,day2,'\');
else
    base_dir = strcat('/cresis/snfs1/data/ICARDS/',date_string(1:4),'/');
    adc_folder_name=strcat(month2,day2,'/');
end
param.year = str2num(date_string(1:4));                                           
param.month = month1;                                      
param.day = day1;                                           
param.radar_name = 'icards';
param.season_name = strcat(date_string(1:4),'_Greenland_P3');                
param.file_regexp = '\S+\.[0-9][0-9][0-9]$';
plot_en = 0; % Set to 1 for gps plots.

if param.year == 1993 || param.year == 1995 || param.year == 1998
  gps_checksum_en = false; % Only 1993,1995,1998? must have this be false
else
  gps_checksum_en = false;
end

% Old ICARDS .mat files have wrong sign for longitude usually, so setting
% this to -1 corrects it
param.longitude_sign = 1;

reuse_tmp_files  = false;

MIN_SEG_SIZE = 1;
MAX_TIME_GAP = 1000/75;
MIN_PRF = 100;

%% Automated Section
full_dir = fullfile(base_dir,adc_folder_name);

fns = get_filenames(full_dir,'','','',struct('regexp',param.file_regexp));
valid_data_file=icards_data_ignore_list(fns,full_dir);%ignore unreasonable raw data files---qishi
fns=fns(valid_data_file);

%If segment goes from 1 to 10 the files for that segment are .000 to .009
for fn_idx = 1:size(fns,1)
  fn = fns{fn_idx};
  [~,fn_name,fn_ext] = fileparts(fn);
  fprintf('Loading %s (%s)\n', fn, datestr(now,'HH:MM:SS'));
  
  tmp_hdr_fn = ct_filename_tmp(param,'','headers',[fn_name fn_ext '.mat']);
  tmp_hdr_fn_dir = fileparts(tmp_hdr_fn);
  if ~exist(tmp_hdr_fn_dir,'dir')
    mkdir(tmp_hdr_fn_dir);
  end
  
  if reuse_tmp_files && exist(tmp_hdr_fn,'file')
    continue;
  end
  
  av = icards_get_available(fn);
  hdr = [];
  
  if any(av(:,1)==10)
    % Corrected Time
    hdr.corr_time = icards_get_data(fn,10);
  end
  
  if any(av(:,1)==11)
    % Corrected Latitude
    hdr.corr_lat = icards_get_data(fn,11);
  end
  
  if any(av(:,1)==12)
    % Corrected Longitude
    hdr.corr_lon = icards_get_data(fn,12);
  end
  
  if any(av(:,1)==13)
    % Corrected Elevation
    hdr.corr_elev = icards_get_data(fn,13);
  end
  
  if any(av(:,1)==5)
    % Computer time
    hdr.comp_time = icards_get_data(fn,5);
  end
  
  if any(av(:,1)==4)
    % NMEA string
    gps_data = char(icards_get_data(fn,4).');
    hdr.nmea_time = NaN*zeros(1,size(gps_data,1));
    hdr.nmea_lat = NaN*zeros(1,size(gps_data,1));
    hdr.nmea_lon = NaN*zeros(1,size(gps_data,1));
    hdr.nmea_elev = NaN*zeros(1,size(gps_data,1));
    hdr.offset = NaN*zeros(1,size(gps_data,1));
    
    for rline = 1:size(gps_data,1)
      if mod(rline-1,1000) == 0
        fprintf('  Range line %d of %d\n', rline, size(gps_data,1));
      end
      C = textscan(gps_data(rline,:),'%s','Delimiter',',');
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
      utc_time_epoch = datenum_to_epoch(datenum(param.year,param.month, ...
        param.day,hour,minute,second));
      hdr.nmea_time(rline) = utc_time_epoch;
      hdr.nmea_lat(rline) = str2double(C{1}{3}(1:2)) + str2double(C{1}{3}(3:end))/60;
      if C{1}{4} ~= 'N'
        hdr.nmea_lat(rline) = -hdr.nmea_lat(rline);
      end
      hdr.nmea_lon(rline) = str2double(C{1}{5}(1:3)) + str2double(C{1}{5}(4:end))/60;
      if C{1}{6} ~= 'E'
        hdr.nmea_lon(rline) = -hdr.nmea_lon(rline);
      end
      hdr.nmea_elev(rline) = str2double(C{1}{10});
      
      hdr.offset(rline)=size(gps_data,2)*rline*2;
   
    end
    
    if all(isnan(hdr.nmea_time))
      warning('only NaN in the file %s,pay attention\n',fn);
    end
    
    save(tmp_hdr_fn,'-struct','hdr');

  end
end

% Load the parsed header data from temporary files
comp_time = [];
corr_time = [];
corr_lat = [];
corr_lon = [];
corr_elev = [];
nmea_time = [];
nmea_lat = [];
nmea_lon = [];
nmea_elev = [];
fn_idxs = [];
for fn_idx = 1:length(fns)
  fn = fns{fn_idx};
  [~,fn_name,fn_ext] = fileparts(fn);
  
  tmp_hdr_fn = ct_filename_tmp(param,'','headers',[fn_name fn_ext '.mat']);
  tmp_hdr_fn_dir = fileparts(tmp_hdr_fn);
  if ~exist(tmp_hdr_fn_dir,'dir')
    mkdir(tmp_hdr_fn_dir);
  end
  
  hdr = load(tmp_hdr_fn);
  
  if isfield(hdr,'nmea_time')
    Nx = length(hdr.nmea_time);
  else
    error('No time field');
  end
  
  if isfield(hdr,'comp_time')
    %For years before 1996, comp_time has different format, use size,1
    %size,2 instead of 1, length
    comp_time =  cat(2,comp_time,reshape(hdr.comp_time,[size(hdr.comp_time,1) size(hdr.comp_time,2)]));
  else
    comp_time = cat(2,comp_time,NaN*ones(1,Nx));
  end
  if isfield(hdr,'corr_time')
    corr_time = cat(2,corr_time,reshape(hdr.corr_time,[1 length(hdr.corr_time)]));
    corr_lat = cat(2,corr_lat,reshape(hdr.corr_lat,[1 length(hdr.corr_lat)]));
    corr_lon = cat(2,corr_lon,reshape(hdr.corr_lon,[1 length(hdr.corr_lon)]));
    corr_elev = cat(2,corr_elev,reshape(hdr.corr_elev,[1 length(hdr.corr_elev)]));
  else
    corr_time = cat(2,corr_time,NaN*ones(1,Nx));
    corr_lat = cat(2,corr_lat,NaN*ones(1,Nx));
    corr_lon = cat(2,corr_lon,NaN*ones(1,Nx));
    corr_elev = cat(2,corr_elev,NaN*ones(1,Nx));
  end
  if isfield(hdr,'nmea_time')
    nmea_time = cat(2,nmea_time,reshape(hdr.nmea_time,[1 length(hdr.nmea_time)]));
    nmea_lat = cat(2,nmea_lat,reshape(hdr.nmea_lat,[1 length(hdr.nmea_lat)]));
    nmea_lon = cat(2,nmea_lon,reshape(hdr.nmea_lon,[1 length(hdr.nmea_lon)]));
    nmea_elev = cat(2,nmea_elev,reshape(hdr.nmea_elev,[1 length(hdr.nmea_elev)]));
  else
    nmea_time = cat(2,nmea_time,NaN*ones(1,Nx));
    nmea_lat = cat(2,nmea_lat,NaN*ones(1,Nx));
    nmea_lon = cat(2,nmea_lon,NaN*ones(1,Nx));
    nmea_elev = cat(2,nmea_elev,NaN*ones(1,Nx));
  end
  
  
  fn_idxs = cat(2,fn_idxs,fn_idx*ones([1 length(hdr.nmea_time)]));
end

%% Fix day wraps in NMEA time
day_wrap_idxs = find(abs(diff(nmea_time) + 86400) < 3600);
for day_wrap_idx = day_wrap_idxs
  warning('Found day wrap, automatically correcting');
  nmea_time(day_wrap_idx+1:end) = nmea_time(day_wrap_idx+1:end) + 86400;
end

%% Get Time Gaps
time_gaps=[];
if all(~isnan(corr_time))
  fprintf('Using corrected time to break segments\n');
  time_gaps = find(abs(diff(corr_time)) > MAX_TIME_GAP);
  time_source = 'corr_time';
elseif all(~isnan(nmea_time))
  fprintf('Using NMEA time to break segments\n');
  time_gaps = find(abs(diff(nmea_time)) > MAX_TIME_GAP);
  time_source = 'nmea_time';
else %if ~all(isnan(corr_time)) && ~all(isnan(nmea_time))
  fprintf('No time source is perfect, handle manually\n');
  % Use prf or pri 
  MAX_TIME_GAP=20;
  all_freq=[1100 2300 4600 9200];
  tmp_nmea_time=nmea_time(~isnan(nmea_time));
  tmp_corr_time=corr_time(~isnan(corr_time));
  %Select which one has more data: corr_time or nmea_time? Or only?
  %if (length(nmea_time)-length(tmp_nmea_time))>(length(corr_time)-length(tmp_corr_time))
  time_source = 'nmea_time';
  tmp_time_gaps=[];
  %get first element NaN
    idx_2=find(~isnan(nmea_time),1);
  for idx_1=idx_2:length(nmea_time)
    if (nmea_time(idx_1)-nmea_time(idx_2)>= MAX_TIME_GAP)
    tmp_time_gaps=[tmp_time_gaps idx_1];
    idx_2=idx_1;
    end   
  end
 
  %Find EPRI
  nmea_gaps=[];
  for tg_idx=1:length(tmp_time_gaps)
    nmea_gaps=[nmea_gaps nmea_time(tmp_time_gaps(tg_idx))];
  end
  
  epri=mean(diff(nmea_gaps))/mean(diff(tmp_time_gaps));%9200/256
  freq=all_freq(4);
  k=1/epri; 
  
  nmea_filled=[];
  %replace nmea NaN values with plausible based on epri
  nmea_ref=find(~isnan(nmea_time),1);
  for idx_fill=1:length(nmea_time)
    if isnan(nmea_time(idx_fill))
    nmea_filled(idx_fill)=nmea_time(nmea_ref)+(idx_fill-nmea_ref)*epri; 
    else
    nmea_ref=idx_fill;
    nmea_filled(idx_fill)=nmea_time(idx_fill);
    end 
  end
  
  %once filled find new time gaps
  nmea_time=nmea_filled;
  time_gaps = find(abs(diff(nmea_time)) > MAX_TIME_GAP);
  time_source = 'nmea_time';
 
%else
 % keyboard;  
  %for this case
end

if plot_en==1
figure(1); clf;
plot(nmea_time - nmea_time(1));
ylabel('GPS time seconds of day');
xlabel('Record');
grid on;
hold on;
plot(time_gaps, nmea_time(time_gaps) - nmea_time(1),'ro');
hold off;
end


%% Break into segments
bad_mask = logical(zeros(size(fns)));
segments = [];
segment_start = fn_idxs(1);
start_time = eval(sprintf('%s(1)', time_source));
seg_idx = 0;
for gap_idx = 1:length(time_gaps)
  time_gap = time_gaps(gap_idx);
  if fn_idxs(time_gap) == fn_idxs(time_gap + 1)
    bad_mask(fn_idxs(time_gap)) = 1;
  end
  
  if bad_mask(fn_idxs(time_gap))
    segment_stop = fn_idxs(time_gap)-1;
  else
    segment_stop = fn_idxs(time_gap);
  end

  if segment_stop - segment_start + 1 >= MIN_SEG_SIZE
    seg_idx = seg_idx + 1;
    segments(seg_idx).start_time = start_time;
    segments(seg_idx).start_idx = segment_start;
    segments(seg_idx).stop_idx = segment_stop;
    [~,fn_start_name,fn_start_name_ext] = fileparts(fns{segment_start});
    [~,fn_stop_name,fn_stop_name_ext] = fileparts(fns{segment_stop});
    fprintf('%2d: %s %3d-%3d %s - %s\n', seg_idx, ...
      datestr(epoch_to_datenum(start_time)), segment_start, segment_stop,...
      [fn_start_name fn_start_name_ext], [fn_stop_name fn_stop_name_ext]);
  end
  
  segment_start = fn_idxs(time_gap)+1;
  start_time = eval(sprintf('%s(time_gap+1)', time_source));
end
segment_stop = length(fns);

if segment_stop - segment_start + 1 >= MIN_SEG_SIZE
  seg_idx = seg_idx + 1;
  segments(seg_idx).start_time = start_time;
  segments(seg_idx).start_idx = segment_start;
  segments(seg_idx).stop_idx = segment_stop;
  [~,fn_start_name,fn_start_name_ext] = fileparts(fns{segment_start});
  [~,fn_stop_name,fn_stop_name_ext] = fileparts(fns{segment_stop});
  fprintf('%2d: %s %3d-%3d %s - %s\n', seg_idx, ...
    datestr(epoch_to_datenum(start_time)), segment_start, segment_stop,...
    [fn_start_name fn_start_name_ext], [fn_stop_name fn_stop_name_ext]);
end
  
[~,sort_idxs] = sort(cell2mat({segments.start_time}));
segments = segments(sort_idxs);

%% Print out some results that can be copied and pasted easily
for seg_idx = 1:length(segments)
  fprintf('%04d%02d%02d\t%02d\t%d\t%d\t%s\t%s\t%s\n', param.year, param.month, param.day, ...
    seg_idx, segments(seg_idx).start_idx, segments(seg_idx).stop_idx, base_dir, adc_folder_name, param.file_regexp);
end
