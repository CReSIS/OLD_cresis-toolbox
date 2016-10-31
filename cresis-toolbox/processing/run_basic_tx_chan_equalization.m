% script run_basic_tx_chan_equalization
%
% Runs basic_tx_chan_equalization
%
% Author: John Paden


% Enable Just One Radar Setup
radar_setup = 'MCORDS5';

%% RDS: MCORDS5
if strcmpi(radar_setup,'MCORDS5')
  [param,defaults] = default_radar_params_2016_Greenland_Polar6_mcords;
  
  % .file_search_mode: Specify how to search for a file: 'last_file',
  %   'specific', 'default', 'segment', 'map', or empty to be asked
  param.file_search_mode = 'segment';

  % .base_dir_search: cell vector of paths to search for data files
  param.base_dir_search = {'D:\awi\','/cresis/snfs1/scratch/2016_Germany_AWI_tests/AWI_ICE_bak/test_flight','/mnt/AWI_SSD0/1604261202/UWB/'};
  
  % out_xml_fn_dir = String containg the directory where the new XML file
  %   will be placed
  if ispc
    param.out_xml_fn_dir = 'C:\waveforms\';
  else
    param.out_xml_fn_dir = '~/waveforms/';
  end

end

%% Automated Section

basic_tx_chan_equalization(param,defaults);
