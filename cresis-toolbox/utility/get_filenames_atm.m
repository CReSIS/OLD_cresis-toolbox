function atm_fns = get_filenames_atm(location,YYYYMMDD,varargin)
% Gets a cell of absolute filename strings for atm L2 smooth nadir files.
%   Only for ATM_smooth_nadir on cresis/projects/metadata/
%
% Input:
%   LOCATION: (string) can be 'arctic' or 'antarctic'
%   YYYYMMDD*: (string) yearmonthday* to get filenames for (ex. '20110331' or '20110331_01')
%   OPTIONAL: (string) data_support_path (default = gRadar.data_support_path)
%
% Output:
%   atm_fns: (cell of strings) atm filenames (same as get_filenames output)
%
% Examples:
%   atm_fns = get_filenames_atm('arctic','20110331')
%   atm_fns = get_filenames_atm('arctic','20110331_01')
%   atm_fns = get_filenames_atm('antarctic','20110331_01_001')
%   with optional data_support_path:
%     atm_fns = get_filenames_atm('arctic','20110331',gRadar.data_support_path)
%     atm_fns = get_filenames_atm('arctic','20110331_01','P:\metadata\')
%     atm_fns = get_filenames_atm('arctic','20110331_01_001','/cresis/projects/metadata/')
%
% Author: Kyle W. Purdon
%
% see also get_filenames.m

% pull out year month day
atmYYYY = YYYYMMDD(1:4);
atmMM = YYYYMMDD(5:6);
atmDD = YYYYMMDD(7:8);

% specify the system correct drive prefix (use data_support_path if possible)
if nargin > 2
  data_support_path = varargin{1};
else
  global gRadar;
  data_support_path = gRadar.data_support_path;
end

% build the base directory and get filenames
if str2double(atmYYYY) >= 2013
  atm_base_dir = fullfile(data_support_path,'ATM_smooth_nadir','ILATM2.002',strcat(atmYYYY,'.',atmMM,'.',atmDD));
  atm_fns = get_filenames(atm_base_dir,'ILATM2','nadir*seg','.csv');
elseif str2double(atmYYYY) >= 2009
  % HANDLE 2009+ ILATM2 ARCTIC/ANTARCTIC CASE
  atm_base_dir = fullfile(data_support_path,'ATM_smooth_nadir','ILATM2.001',strcat(atmYYYY,'.',atmMM,'.',atmDD));
  atm_fns = get_filenames(atm_base_dir,'ILATM2','nadir*seg','');
elseif str2double(atmYYYY) <= 2008 && strcmp(location,'antarctic')
  % HANDLE 1993-2008 ANTARCTIC CASE
  atm_base_dir = fullfile(data_support_path,'ATM_smooth_nadir','BLATM2_ATMicessn_v01',strcat(atmYYYY,'_AN_NASA'),YYYYMMDD(1:8));
  atm_fns = get_filenames(atm_base_dir,'','nadir*seg','');
elseif str2double(atmYYYY) <= 2007 && strcmp(location,'arctic')
  % HANDLE 1993-2007 ARCTIC CASE
  directory_add_on = {'','a','b','c','_4cT3'};
  atm_fns = {};
  for idx = 1:4
    atm_base_dir = fullfile(data_support_path,'ATM_smooth_nadir','BLATM2_ATMicessn_v01','icessn',strcat('AIM',atmYYYY),strcat(YYYYMMDD(1:8),directory_add_on{idx}));
    atm_fns = cat(1,atm_fns,get_filenames(atm_base_dir,'','nadir*seg',''));
  end
elseif str2double(atmYYYY) == 2008 && strcmp(location,'arctic')
  % HANDLE 2008 ARCTIC CASE
  atm_base_dir = fullfile(data_support_path,'ATM_smooth_nadir','BLATM2_ATMicessn_v01','icessn',strcat('AIM',atmYYYY),strcat('ATM_',YYYYMMDD(3:8)));
  atm_fns = get_filenames(atm_base_dir,'','nadir*seg','');
else
  error('YYYYMMDD should be a string, no year could be determined from given input.');
end

% WARN IF ATM_FNS IS EMPTY
if isempty(atm_fns)
  warning('No ATM filenames found, atm_fns is empty');
end

end