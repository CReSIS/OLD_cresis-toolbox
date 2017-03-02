% script run_basic_noise_analysis
%
% Runs basic_noise_analysis
%
% Author: John Paden

% Initialization
physical_constants;
param = [];

% Enable Just One Radar Setup
radar_setup = 'MCORDS3_DC8';

%% RDS: MCORDS5
if strcmpi(radar_setup,'MCORDS5')
  [param,defaults] = default_radar_params_2016_Greenland_Polar6_mcords;
  
  % .file_search_mode: Specify how to search for a file: 'last_file',
  %   'specific', 'default', or empty to be asked
  param.file_search_mode = 'last_file';

  % .base_dir_search: cell vector of paths to search for data files
  param.base_dir_search = {'D:\awi\','/mnt/AWI_SSD0/1604261101/UWB/','/mnt/AWI_SSD0/1604261202/UWB/'};
  
  % .pdf_en: Enable time domain, probability density function, and quantization plots
  param.pdf_en = false;
  % .psd_en: Enable PSD plot
  param.psd_en = true;
  
  % .img: wf-adc pair list which specifies which waveform-adc pairs to
  %   analyze
  wf = 2; adcs = 1:24;
  param.img = cat(2, wf*ones(24,1), adcs.');

  % .recs: two element vector specifying which records/range-lines to load
  %   [start_record num_records]
  param.recs = [0 250];
  
  % .presums: Number of additional software presums (coherent averaging) to do
  param.presums = 1;
end

%% RDS: MCORDS3_DC8
if strcmpi(radar_setup,'MCORDS3_DC8')
  [param,defaults] = default_radar_params_2016_Antarctica_DC8_mcords;
  
  % .file_search_mode: Specify how to search for a file: 'last_file',
  %   'specific', 'default', or empty to be asked
  param.file_search_mode = 'segment';

  % .base_dir_search: cell vector of paths to search for data files
  param.base_dir_search = {'D:\awi\','/process/mcords/','/mnt/AWI_SSD0/1604261202/UWB/'};
  
  % .pdf_en: Enable time domain, probability density function, and quantization plots
  param.pdf_en = false;
  % .psd_en: Enable PSD plot
  param.psd_en = true;
  
  % .img: wf-adc pair list which specifies which waveform-adc pairs to
  %   analyze
  wf = 3; adcs = 1:6;
  param.img = cat(2, wf*ones(6,1), adcs.');

  % .recs: two element vector specifying which records/range-lines to load
  %   [start_record num_records]
  param.recs = [];
  
  % .presums: Number of additional software presums (coherent averaging) to do
  param.presums = 1;
end

%% Automated Section

basic_noise_analysis(param,defaults);

return;
