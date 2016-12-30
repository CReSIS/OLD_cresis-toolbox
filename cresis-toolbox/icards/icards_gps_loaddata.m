%NOTES:This function load NMEA and TRAJ data of one certain day !  
today_string=num2str(today);
param1.date=num2str(today_string(1:8));
param1.year=num2str(today_string(3:4));
param.year = str2num(today_string(1:4)); 
if ispc
  out_fn='X:\metadata\';
  out_dir=fullfile(out_fn,which_season(11:27),'\');
  base_dir=strcat('Z:\ICARDS\',num2str(param.year),'\');
  in_base_path = 'Z:\NASA\';
else
  out_fn='/cresis/snfs1/dataproducts/metadata/';
  out_dir=fullfile(out_fn,which_season(11:27),'/');
  base_dir=strcat('/cresis/snfs1/data/ICARDS/',num2str(param.year),'/');
  in_base_path = '/cresis/snfs1/data/NASA';
end
  
  full_file_out=strcat(out_dir,param1.date,'_nmea','.csv');
  D= datevec(param1.date,'yyyymmdd');
  day= daysact(strcat('1-1-',num2str(param1.year)),sprintf('%4d-%0.2d-%2d',D(1),D(2),D(3)));
  
  [month1,month2,day1,~]=icards_monthANDday(today);
  
  adc_folder_name=strcat(month2,today_string(7:8)); 
  param.month = month1;
  param.day = day1;                                   
  param.radar_name = 'icards';                    
  param.season_name=strcat(num2str(param.year),'_Greenland_P3');
  param.file_regexp = '\S+\.[0-9][0-9][0-9]$';

  plot_en = 0; % Set to 1 for gps plots.

  if param.year == 1993 || param.year == 1995
    gps_checksum_en = false; % Only 1993,1995 must have this be false
  else
    gps_checksum_en = true; % Only 1993,1995 must have this be false
  end

  % Old ICARDS .mat files have wrong sign for longitude usually, so setting
  % this to -1 corrects it
  param.longitude_sign = 1;

  reuse_tmp_files  = true;

  MIN_SEG_SIZE = 1;
  MAX_TIME_GAP = 1000/75;
  MIN_PRF = 100;

  %%Automated Section
  full_dir = fullfile(base_dir,adc_folder_name);
  fns = get_filenames(full_dir,'','','',struct('regexp',param.file_regexp));
  valid_data_file=icards_data_ignore_list(fns,full_dir);%ignore unreasonable temporary files---qishi
  fns=fns(valid_data_file);
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
    hdr.nmea_time = hdr.nmea_time;
    hdr.nmea_lat = hdr.nmea_lat;
    hdr.nmea_lon = hdr.nmea_lon;
    hdr.nmea_elev = hdr.nmea_elev;

    if (all(isnan(hdr.nmea_time)))||(all(isnan(hdr.nmea_lat)))||(all(isnan(hdr.nmea_lon)))||(all(isnan(hdr.nmea_elev)))%to ingnore this temporary file which contains only NaN---qishi
      warning('no valid data in %s except NaN,to load next file if there are any\n',fn);
      continue;
    end

    if isfield(hdr,'nmea_time')
      Nx = length(hdr.nmea_time);
    else
      error('No time field');
    end

    if isfield(hdr,'comp_time')
      comp_time = cat(2,comp_time,reshape(hdr.comp_time(1,:),[1 length(hdr.comp_time(1,:))]));
    else
      comp_time = cat(2,comp_time,NaN*ones(1,Nx));
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


    if (all(isnan(nmea_time)))||(all(isnan(nmea_lat)))||(all(isnan(nmea_lon)))||(all(isnan(nmea_elev)));%to ingnore this temporary file which contains only NaN---qishi
      warning('No valid data in %s except NaN,to load next file if there are any',fn);%---qishi
      continue;
    end
  end

  cur_time = nmea_time(1);
  bad_mask = logical(zeros(size(nmea_time)));
  if any(diff(nmea_time)<=0)
    warning('Repeated time or decreased time point exists in NMEA time,jump those points');%---qishi
  end
  for idx = 2:length(nmea_time)
    if nmea_time(idx) <= cur_time
      bad_mask(idx) = 1;
    else
      cur_time = nmea_time(idx);
    end
  end
  nmea_time = nmea_time(~bad_mask);
  nmea_elev = nmea_elev(~bad_mask);
  nmea_lat = nmea_lat(~bad_mask);
  nmea_lon = nmea_lon(~bad_mask);

  NMEA.time=nmea_time+utc_leap_seconds(nmea_time(1));%add leap seconds here to get "gps time"---qishi
  NMEA.elev=nmea_elev;
  NMEA.lat=nmea_lat;
  NMEA.lon=nmea_lon;
  NMEA.roll=zeros(1,length(NMEA.time));
  NMEA.pitch=zeros(1,length(NMEA.time));
  NMEA.heading=zeros(1,length(NMEA.time));

  if any(isnan(NMEA.time))%there may be NaNs in temporary files---qishi
    warning('NaN exists in NMEA time,jump those points\n');
  end
  nan_mask=isnan(NMEA.time);
  NMEA.time=NMEA.time(~nan_mask);
  NMEA.elev=NMEA.elev(~nan_mask);
  NMEA.lat=NMEA.lat(~nan_mask);
  NMEA.lon=NMEA.lon(~nan_mask);
  NMEA.roll=NMEA.roll(~nan_mask);
  NMEA.pitch=NMEA.pitch(~nan_mask);
  NMEA.heading=NMEA.heading(~nan_mask);

  year=str2num(today_string(1:4));
  if year<1999
    traj_mark_in=0;
  else
    traj_mark_in=1;
  end
  [traj_mark_out,season_name,traj_folder,traj_filename,gps_source]=icards_data_traj_list(strcat(today_string,'_01'),traj_mark_in);% for the convenience to load traj files quickly---qishi
  if traj_mark_out % not all days of a certain seanson has trajactory file---qishi      
    in_fn = fullfile(in_base_path,traj_folder,traj_filename);
    file_type = 'Traj';
    params = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc');
    sync_flag = 0;

    gps_tmp = read_gps_traj(in_fn,params);
    idx_not_nan=find(isnan(gps_tmp.gps_time)==0);
    gps_tmp.lat=gps_tmp.lat(idx_not_nan);
    gps_tmp.lon=gps_tmp.lon(idx_not_nan);
    gps_tmp.elev=gps_tmp.elev(idx_not_nan);
    gps_tmp.gps_time=gps_tmp.gps_time(idx_not_nan);
    gps_tmp.roll=gps_tmp.roll(idx_not_nan);
    gps_tmp.pitch=gps_tmp.pitch(idx_not_nan);
    gps_tmp.heading=gps_tmp.heading(idx_not_nan);
    TRAJ=gps_tmp;%TRAJ Data
    fprintf('Trajactory file is available for %d,loading file %s\n',today,in_fn);
  else
    TRAJ=[];
  end

