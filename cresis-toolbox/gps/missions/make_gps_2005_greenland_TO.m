% script make_gps_2005_greenland_TO_GPS
%
% Makes the GPS files for 2005 Greenland Twin Otter field season

tic;

global gRadar;

support_path = '';
data_support_path = '/cresis/snfs1/dataproducts/metadata/lsmith/';

if isempty(support_path)
  support_path = gRadar.support_path;
end

gps_path = fullfile(support_path,'gps','2005_Greenland_TO');
if ~exist(gps_path,'dir')
  fprintf('Making directory %s\n', gps_path);
  fprintf('  Press a key to proceed\n');
  pause;
  mkdir(gps_path);
end

if isempty(data_support_path)
  data_support_path = gRadar.data_support_path;
end

% ======================================================================
% User Settings
% ======================================================================
debug_level = 1;

in_base_path = fullfile(data_support_path,'NMEA/2005_Greenland_TO');
% in_base_path = fullfile(data_support_path);

file_idx = 0; in_fns = {}; out_fns = {}; file_type = {}; params = {};
gps_source_to_use = 'NMEA';
if strcmp(gps_source_to_use,'Novatel')
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'rover_LC_ppp_Greenland_20110324.gps');
%   out_fns{file_idx} = 'gps_20110324.mat';
%   file_type{file_idx} = 'Novatel';
%   params{file_idx} = struct('year',2011,'month',3,'day',24,'time_reference','gps');
%   gps_source{file_idx} = 'Novatel-field_20110324-ppp';
  
  
elseif strcmp(gps_source_to_use,'NMEA')
%   
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'');
%   out_fns{file_idx} = 'gps_20050421.mat';
%   file_type{file_idx} = 'NMEA';
%   params{file_idx} = struct('year',2005,'month',04,'day',21,'time_reference','utc','nmea_tag','$GPGGA','format',2,'combine',1);
%   gps_source{file_idx} = 'NMEA';
%   
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'');
%   out_fns{file_idx} = 'gps_20050504.mat';
%   file_type{file_idx} = 'NMEA';
%   params{file_idx} = struct('year',2005,'month',05,'day',04,'time_reference','utc','nmea_tag','$GPGGA','format',2,'combine',1);
%   gps_source{file_idx} = 'NMEA-field';
%   
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'');
%   out_fns{file_idx} = 'gps_20050505.mat';
%   file_type{file_idx} = 'NMEA';
%   params{file_idx} = struct('year',2005,'month',05,'day',05,'time_reference','utc','nmea_tag','$GPGGA','format',2,'combine',1);
%   gps_source{file_idx} = 'NMEA-field';
  
  file_idx = file_idx + 1;
  in_fns{file_idx} = fullfile(in_base_path,'');
  out_fns{file_idx} = 'gps_20050506.mat';
  file_type{file_idx} = 'NMEA';
  params{file_idx} = struct('year',2005,'month',05,'day',06,'time_reference','utc','nmea_tag','$GPGGA','format',2,'combine',1);
  gps_source{file_idx} = 'NMEA-field';
  
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'');
%   out_fns{file_idx} = 'gps_20050507.mat';
%   file_type{file_idx} = 'NMEA';
%   params{file_idx} = struct('year',2005,'month',05,'day',07,'time_reference','utc','nmea_tag','$GPGGA','format',2,'combine',1);
%   gps_source{file_idx} = 'NMEA-field';
%   
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'');
%   out_fns{file_idx} = 'gps_20050508.mat';
%   file_type{file_idx} = 'NMEA';
%   params{file_idx} = struct('year',2005,'month',05,'day',08,'time_reference','utc','nmea_tag','$GPGGA','format',2,'combine',1);
%   gps_source{file_idx} = 'NMEA-field';
%   
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'');
%   out_fns{file_idx} = 'gps_20050511.mat';
%   file_type{file_idx} = 'NMEA';
%   params{file_idx} = struct('year',2005,'month',05,'day',11,'time_reference','utc','nmea_tag','$GPGGA','format',2,'combine',1);
%   gps_source{file_idx} = 'NMEA-field';
%   
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'');
%   out_fns{file_idx} = 'gps_20050512.mat';
%   file_type{file_idx} = 'NMEA';
%   params{file_idx} = struct('year',2005,'month',05,'day',12,'time_reference','utc','nmea_tag','$GPGGA','format',2,'combine',1);
%   gps_source{file_idx} = 'NMEA-field';
%   
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'');
%   out_fns{file_idx} = 'gps_20050513.mat';
%   file_type{file_idx} = 'NMEA';
%   params{file_idx} = struct('year',2005,'month',05,'day',13,'time_reference','utc','nmea_tag','$GPGGA','format',2,'combine',1);
%   gps_source{file_idx} = 'NMEA-field';
%   
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'');
%   out_fns{file_idx} = 'gps_20050514.mat';
%   file_type{file_idx} = 'NMEA';
%   params{file_idx} = struct('year',2005,'month',05,'day',14,'time_reference','utc','nmea_tag','$GPGGA','format',2,'combine',1);
%   gps_source{file_idx} = 'NMEA-field';
%   
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'');
%   out_fns{file_idx} = 'gps_20050517.mat';
%   file_type{file_idx} = 'NMEA';
%   params{file_idx} = struct('year',2005,'month',05,'day',17,'time_reference','utc','nmea_tag','$GPGGA','format',2,'combine',1);
%   gps_source{file_idx} = 'NMEA-field';
%   
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'');
%   out_fns{file_idx} = 'gps_20050518.mat';
%   file_type{file_idx} = 'NMEA';
%   params{file_idx} = struct('year',2005,'month',05,'day',18,'time_reference','utc','nmea_tag','$GPGGA','format',2,'combine',1);
%   gps_source{file_idx} = 'NMEA-field';
%   
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'');
%   out_fns{file_idx} = 'gps_20050519.mat';
%   file_type{file_idx} = 'NMEA';
%   params{file_idx} = struct('year',2005,'month',05,'day',19,'time_reference','utc','nmea_tag','$GPGGA','format',2,'combine',1);
%   gps_source{file_idx} = 'NMEA-field';
%   
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'');
%   out_fns{file_idx} = 'gps_20050522.mat';
%   file_type{file_idx} = 'NMEA';
%   params{file_idx} = struct('year',2005,'month',05,'day',22,'time_reference','utc','nmea_tag','$GPGGA','format',2,'combine',1);
%   gps_source{file_idx} = 'NMEA-field';
%   
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'');
%   out_fns{file_idx} = 'gps_20050523.mat';
%   file_type{file_idx} = 'NMEA';
%   params{file_idx} = struct('year',2005,'month',05,'day',23,'time_reference','utc','nmea_tag','$GPGGA','format',2,'combine',1);
%   gps_source{file_idx} = 'NMEA-field';
  
 
end

% ======================================================================
% Read and translate files according to user settings
% ======================================================================
make_gps;

