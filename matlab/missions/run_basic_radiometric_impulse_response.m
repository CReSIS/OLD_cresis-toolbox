% script run_basic_radiometric_impulse_response
%
% Runs basic_radiometric_impulse_response
%
% Author: John Paden

% Initialization
physical_constants;
param = [];

% Enable Just One Radar Setup
radar_setup = 'MCORDS5';

%% RDS: MCORDS5
if strcmpi(radar_setup,'MCORDS5')
  [param,defaults] = default_radar_params_2016_Greenland_Polar6_mcords;
  
  % .file_search_mode: Specify how to search for a file: 'last_file',
  %   'specific', 'default', or empty to be asked
  param.file_search_mode = '';

  % .base_dir_search: cell vector of paths to search for data files
  param.base_dir_search = {'D:\awi\','/mnt/AWI_SSD0/1604261101/UWB/','/mnt/AWI_SSD0/1604261202/UWB/'};
  
  % .img: wf-adc pair list which specifies which waveform-adc pairs to
  %   analyze
%   param.img = cat(2,1*ones(24,1),[1:24].'); % First waveform
%   param.img = cat(2,2*ones(24,1),[1:24].'); % Second waveform
  param.img = cat(2,3*ones(24,1),[1:24].'); % Third waveform

  % .recs: two element vector specifying which records/range-lines to load
  %   [start_record num_records]
  param.recs = [0 inf];
  
  % .presums: Number of additional software presums (coherent averaging) to do
  param.presums = 1;
end

%% RDS: MCORDS3_P3
if strcmpi(radar_setup,'MCORDS3_P3')
  [param,defaults] = default_radar_params_2014_Greenland_P3_mcords;
  
  % .file_search_mode: Specify how to search for a file: 'last_file', 'specific'
  param.file_search_mode = '';

  % .base_dir_search: cell vector of paths to search for data files
  param.base_dir_search = {'/cresis/snfs1/data/MCoRDS/2014_Greenland_P3/20140401/'};
  
  % .img: wf-adc pair list which specifies which waveform-adc pairs to
  %   analyze
  param.img = cat(2,2*ones(15,1),[2:16].'); param.ref_wf_adc = 3;
  
  % .recs: two element vector specifying which records/range-lines to load
  %   [start_record num_records]
  param.recs = [0 inf];
  
  % .presums: Number of additional software presums (coherent averaging) to do
  param.presums = 1;
end



%% Automated Section

basic_radiometric_impulse_response(param,defaults);

return;
