% script make_gps_2015_greenland_LC130
%
% Makes the GPS files for 2015 Greenland LC130 field season

tic;

global gRadar;

support_path = '';
data_support_path = '';

if isempty(support_path)
  support_path = gRadar.support_path;
end

gps_path = fullfile(support_path,'gps','2015_Greenland_LC130');
if ~exist(gps_path,'dir')
  fprintf('Making directory %s\n', gps_path);
  fprintf('  Press a key to proceed\n');
  pause;
  mkdir(gps_path);
end

if isempty(data_support_path)
  data_support_path = gRadar.data_support_path;
end

debug_level = 1;

in_base_path = fullfile(data_support_path,'2015_Greenland_LC130');

file_idx = 0; in_fns = {}; out_fns = {}; file_type = {}; params = {}; gps_source = {};
sync_fns = {}; sync_params = {};

% gps_source_to_use = 'NMEA';
gps_source_to_use = 'sonntag';
% gps_source_to_use = 'ATM-field';
% gps_source_to_use = 'ATM';

if strcmpi(gps_source_to_use,'NMEA')
%   file_idx = file_idx + 1;
%   year = 2015; month = 3; day = 24;
%   in_fns{file_idx} = get_filenames(in_base_path,sprintf('GPS_%04d%02d%02d',year,month,day),'','.txt');
%   out_fns{file_idx} = sprintf('gps_%04d%02d%02d.mat',year,month,day);
%   file_type{file_idx} = 'NMEA';
%   params{file_idx} = struct('year',year,'month',month,'day',day,'format',1,'time_reference','utc');
%   gps_source{file_idx} = 'nmea-field';
%   sync_flag{file_idx} = 0;

  file_idx = file_idx + 1;
  year = 2015; month = 3; day = 25;
  in_fns{file_idx} = get_filenames(in_base_path,sprintf('GPS_%04d%02d%02d',year,month,day),'','.txt');
  out_fns{file_idx} = sprintf('gps_%04d%02d%02d.mat',year,month,day);
  file_type{file_idx} = 'NMEA';
  params{file_idx} = struct('year',year,'month',month,'day',day,'format',1,'time_reference','utc');
  gps_source{file_idx} = 'nmea-field';
  sync_flag{file_idx} = 0;
  
elseif strcmpi(gps_source_to_use,'sonntag')

  file_idx = file_idx + 1;
  year = 2015; month = 3; day = 25;
  in_fns{file_idx} = get_filenames(in_base_path,sprintf('%02d%02d%02d',mod(year,100),month,day),'','.gps');
  out_fns{file_idx} = sprintf('gps_%04d%02d%02d.mat',year,month,day);
  file_type{file_idx} = 'traj';
  params{file_idx} = struct('headerlines',0,'year',year,'time_reference','gps');
  gps_source{file_idx} = 'sonntag-field';
  sync_flag{file_idx} = 0;
  
elseif strcmpi(gps_source_to_use,'ATM-field')
%   file_idx = file_idx + 1;
%   year = 2015; month = 3; day = 13;
%   in_fns{file_idx} = fullfile(in_base_path,'BD982_13Mar15_PPPK_P15Mar15_FROM_DMS.out');
%   out_fns{file_idx} = sprintf('gps_%04d%02d%02d.mat',year,month,day);
%   file_type{file_idx} = 'applanix';
%   params{file_idx} = struct('year',year,'month',month,'day',day,'format',3,'time_reference','gps');
%   gps_source{file_idx} = 'atm-field';
%   sync_flag{file_idx} = 0;
  
%   file_idx = file_idx + 1;
%   year = 2015; month = 3; day = 19;
%   in_fns{file_idx} = fullfile(in_base_path,'BD982_19Mar15_GNSSK_P19Mar15.out');
%   out_fns{file_idx} = sprintf('gps_%04d%02d%02d.mat',year,month,day);
%   file_type{file_idx} = 'applanix';
%   params{file_idx} = struct('year',year,'month',month,'day',day,'format',3,'time_reference','utc');
%   gps_source{file_idx} = 'atm-field';
%   sync_flag{file_idx} = 0;
   
  file_idx = file_idx + 1;
  year = 2015; month = 3; day = 24;
  in_fns{file_idx} = fullfile(in_base_path,'BD960_24Mar15_GNSSK_P24Mar15_V5eth.out');
  out_fns{file_idx} = sprintf('gps_%04d%02d%02d.mat',year,month,day);
  file_type{file_idx} = 'applanix';
  params{file_idx} = struct('year',year,'month',month,'day',day,'format',3,'time_reference','utc');
  gps_source{file_idx} = 'atm-field';
  sync_flag{file_idx} = 0; 
  
elseif strcmpi(gps_source_to_use,'ATM')
  % Just some simple code to automate creation of the code in this section:
  %
  % ATM_fns = get_filenames(in_base_path,'','','.out');
  % fn_dates = [];
  % for idx = 1:length(ATM_fns)
  %   fn = ATM_fns{idx};
  %   [~,fn_name] = fileparts(fn);
  %   fn_dates(end+1) = datenum(sprintf('%s %s, 20%s', fn_name(9:11), fn_name(7:8), fn_name(12:13)));
  % end
  % fn_dates = sort(fn_dates);
  % for idx = 1:length(fn_dates)
  %   [year,month,day] = datevec(fn_dates(idx));
  %   fprintf('year = %d; month = %d; day = %d;\n', year, month, day);
  % end
  
end


% ======================================================================
% Read and translate files according to user settings
% ======================================================================
make_gps;




