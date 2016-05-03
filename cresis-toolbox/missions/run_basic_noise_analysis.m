% script run_basic_noise_analysis
%
% Runs basic_noise_analysis
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
  param.file_search_mode = 'last_file';

  % .base_dir_search: cell vector of paths to search for data files
  param.base_dir_search = {'D:\awi\chan1','/mnt/HDD6/1604130301/UWB/chan1/','/mnt/AWI_SSD0/1604180702/UWB/chan1/'};
  
  % .pdf_en: Enable time domain, probability density function, and quantization plots
  param.pdf_en = false;
  % .psd_en: Enable PSD plot
  param.psd_en = true;
  
  % .noise_rbins: Specify which bins to use for noise analysis (nonpositive
  %   numbers cause basic_noise_analysis to reference from the end of the range line)
  param.noise_rbins = [-499 0];
  
  % .img: wf-adc pair list which specifies which waveform-adc pairs to
  %   analyze
  param.img = cat(2,3*ones(24,1),[1:24].');
  % param.img = cat(2,3*ones(8,1),[9:16].'); % Center subarray
%   param.img = cat(2,3*ones(8,1),[1:8].'); % Left subarray
  % param.img = cat(2,3*ones(8,1),[17:24].'); % Right subarray

  % .recs: two element vector specifying which records/range-lines to load
  %   [start_record num_records]
  param.recs = [0 250];
  
  % .presums: Number of additional software presums (coherent averaging) to do
  param.presums = 1;
end

%% Automated Section

basic_noise_analysis(param,defaults);

return;
