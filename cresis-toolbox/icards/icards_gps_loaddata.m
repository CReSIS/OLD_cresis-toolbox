%NOTES:This function load NMEA and TRAJ data of one certain day !
%% Preparations
param1 = read_param_xls(ct_filename_param('rds_param_2002_Greenland_P3.xls'),'20020520_01');%change the data (segment does not matter)
param1.date=param1.day_seg(1:8);
param1.year=param1.day_seg(3:4);
out_fn='X:\metadata\';
out_dir=fullfile(out_fn,param1.season_name,'\');
full_file_out=strcat(out_dir,param1.date,'_nmea','.csv');
D= datevec(param1.date,'yyyymmdd');
day= daysact(strcat('1-1-',num2str(param1.year)),sprintf('%4d-%0.2d-%2d',D(1),D(2),D(3)));
base_dir = 'Z:\ICARDS\2002\';                                %change needed
adc_folder_name='may20\';                                    %change needed
param.year = 2002;                                           %change needed
param.month = 5;
param.day = 20;                                              %change needed
param.radar_name = 'icards';
param.season_name = '2002_Greenland_P3';                     %change needed
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
% fns = fns(2:205);  %to ignore "fiberdly" of 20020524!!!!!!!!!!!!!!!!!!!!!
fns = fns(1:43);   %to ignore "seaice"   of 20020520!!!!!!!!!!!!!!!!!!!!!

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
   
  if isfield(hdr,'nmea_time')
    Nx = length(hdr.nmea_time);
  else
    error('No time field');
  end
  
  if isfield(hdr,'comp_time')
    comp_time = cat(2,comp_time,reshape(hdr.comp_time,[1 length(hdr.comp_time)]));
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
end

cur_time = nmea_time(1);
bad_mask = logical(zeros(size(nmea_time)));
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

NMEA.time=nmea_time;
NMEA.elev=nmea_elev;
NMEA.lat=nmea_lat;
NMEA.lon=nmea_lon;
NMEA.roll=zeros(1,length(NMEA.time));
NMEA.pitch=zeros(1,length(NMEA.time));
NMEA.heading=zeros(1,length(NMEA.time));

if 0                                     %ONLY enable when traj file exists
%% Load TRAJ FILES 
global gRadar;
data_support_path = '';

if isempty(data_support_path)
  data_support_path = gRadar.data_support_path;
end

debug_level = 1;

in_base_path = fullfile(data_support_path,'2002_Greenland_P3');%change needed

in_fn = fullfile(in_base_path,'2002_Greenland_P3_traj','020524_aa_l12_jgs_itrf00_17jul02_npm');%change needed
file_type = 'Traj';
params = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc');
gps_source = 'atm-final_20020524';                           %change needed
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
end

