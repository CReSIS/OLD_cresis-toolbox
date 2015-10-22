% script make_gps_2002_greenland_P3
% Makes the DGPSwINS????? files for 2002 Greenland P3 field season
%see icards_gps_missinNASA_csv.m to get csv files for days without
%trajectory data. (check time reference: should be gps)
tic;

global gRadar;

support_path = '';
data_support_path = '';

if isempty(support_path)
  support_path = gRadar.support_path;
end

gps_path = fullfile(support_path,'gps','2002_Greenland_P3');
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

in_base_path = fullfile(data_support_path,'2002_Greenland_P3');
gps_path = fullfile(support_path,'gps','2002_Greenland_P3');

file_idx = 0; in_fns = {}; out_fns = {}; file_type = {}; params = {};

% file_idx = file_idx + 1;
% in_fns{file_idx} = fullfile(in_base_path,'sbet_27MarPPPFinal.out');
% out_fns{file_idx} = 'gps_20090327.mat';
% file_type{file_idx} = 'Applanix';
% params{file_idx} = struct('year',2009,'month',03,'day',27,'time_reference','utc');
% gps_source{file_idx} = 'atm-final_20090717';

% file_idx = file_idx + 1;
% in_fns{file_idx} = fullfile(in_base_path,'sbet_30MarPPPFinal.out');
% out_fns{file_idx} = 'gps_20090330.mat';
% file_type{file_idx} = 'Applanix';
% params{file_idx} = struct('year',2009,'month',03,'day',30,'time_reference','utc');
% gps_source{file_idx} = 'atm-final_20090717';

%%%%%%%%%%%%%%%%%%
%2001 Greenland_P3

%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'05-19-1_gps.csv');
%   out_fns{file_idx} = 'gps_20010519.mat';
%   file_type{file_idx} = 'csv';
%   params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','gps');
%   gps_source{file_idx} = 'atm-final_20010519';
%   
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'05-20-1_gps.csv');
%   out_fns{file_idx} = 'gps_20010520.mat';
%   file_type{file_idx} = 'csv';
%   params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','gps');
%   gps_source{file_idx} = 'atm-final_20010520';
  
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'2001_Greenland_P3_traj','010521_aa_l12_jgs_itrf97_19jul01');
%   out_fns{file_idx} = 'gps_20010521.mat';
%   file_type{file_idx} = 'Traj';
%   params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc');
%   gps_source{file_idx} = 'atm-final_20010521';
%   sync_flag{file_idx} = 0;
%   
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'2001_Greenland_P3_traj','010523_aa_l12_jgs_itrf97_09aug01');
%   out_fns{file_idx} = 'gps_20010523.mat';
%   file_type{file_idx} = 'Traj';
%   params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc');
%   gps_source{file_idx} = 'atm-final_20010523';
%   
%   file_idx = file_idx + 1;
%   in_fns{file_idx} = fullfile(in_base_path,'2001_Greenland_P3_traj','010524_aa_l12_jgs_itrf97_06jul01');
%   out_fns{file_idx} = 'gps_20010524.mat';
%   file_type{file_idx} = 'Traj';
%   params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc');
%   gps_source{file_idx} = 'atm-final_20010524';
  
  file_idx = file_idx + 1;
  in_fns{file_idx} = fullfile(in_base_path,'20020604_nmea.csv');
  out_fns{file_idx} = 'gps_20020604.mat';
  file_type{file_idx} = 'csv';
  params{file_idx} = struct('input_format','%f%f%f%f%f%f%f%f%f%f%f%f%f','time_reference','utc','type',[3]);%add a new type valued 3 for "read_gps_csv" to process
  gps_source{file_idx} = 'atm-final_20020604'; 
 %%%%%%%%%%%%%%%%%%%%

% file_idx = file_idx + 1;
% in_fns{file_idx} = fullfile(in_base_path,'sbet_28AprPPPFinal.out');
% out_fns{file_idx} = 'gps_20090428.mat';
% file_type{file_idx} = 'Applanix';
% params{file_idx} = struct('year',2009,'month',04,'day',28,'time_reference','utc');
% gps_source{file_idx} = 'atm-final_20090717';
% 
% file_idx = file_idx + 1;
% in_fns{file_idx} = fullfile(in_base_path,'sbet_01MayPPPFinal.out');
% out_fns{file_idx} = 'gps_20090501.mat';
% file_type{file_idx} = 'Applanix';
% params{file_idx} = struct('year',2009,'month',05,'day',01,'time_reference','utc');
% gps_source{file_idx} = 'atm-final_20090717';
% 
% file_idx = file_idx + 1;
% in_fns{file_idx} = fullfile(in_base_path,'sbet_02MayPPPFinal.out');
% out_fns{file_idx} = 'gps_20090502.mat';
% file_type{file_idx} = 'Applanix';
% params{file_idx} = struct('year',2009,'month',05,'day',01,'time_reference','utc');
% gps_source{file_idx} = 'atm-final_20090717';
% 
% file_idx = file_idx + 1;
% in_fns{file_idx} = fullfile(in_base_path,'sbet_06MayFlt1_PPPFinal.out');
% out_fns{file_idx} = 'gps_20090506A.mat';
% file_type{file_idx} = 'Applanix';
% params{file_idx} = struct('year',2009,'month',05,'day',06,'time_reference','utc');
% gps_source{file_idx} = 'atm-final_20090717';
% 
% file_idx = file_idx + 1;
% in_fns{file_idx} = fullfile(in_base_path,'sbet_06MayFlt2_PPPFinal.out');
% out_fns{file_idx} = 'gps_20090506B.mat';
% file_type{file_idx} = 'Applanix';
% params{file_idx} = struct('year',2009,'month',05,'day',06,'time_reference','utc');
% gps_source{file_idx} = 'atm-final_20090717';

% ======================================================================
% Read and translate files according to user settings
% ======================================================================
make_gps;

