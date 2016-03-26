% script cresis_xml_mapping
%
% XML variable name mapping for each version of the NI XML files. The NI
% XML files contain the radar settings. They do not always correspond
% to the actual settings used by the radar due to bugs, but are correct
% most of the time. Version 2.0 onward files should not have these bugs.
%
% https://wiki.cresis.ku.edu/cresis/National_Instruments_XML_File_Guide
%
% Inputs:
%  xml_version: double scalar. Valid options are 1.0, 1.1, 1.2, 1.3, 2.0
%
% Outputs:
%  xml_file_prefix (will be created if it does not exist)
%
% Author: John Paden

if xml_version == 1.0
  config_var = 'Configuration';
  config_var_enc = 'Configuration';
  prf_var = 'PRF';
  prf_var_enc = 'PRF';
  ram_var = 'RAM_Taper';
  ram_var_enc = 'RAMZ20Taper';
  ram_amp_var = 'Ram_Amplitude';
  ram_amp_var_enc = 'RamZ20Amplitude';
  phase_var = 'Phase_Offset_deg';
  phase_var_enc = 'PhaseZ20Offset';
  wave_var_enc = 'Z23Wave';
  ttl_start_var_enc = 'TTLZ20Start';
  ttl_length_var_enc = 'TTLZ20Length';
  if ~exist('xml_file_prefix','var')
    xml_file_prefix = 'DDS';
  end
elseif xml_version == 1.1
  config_var = 'Configuration';
  config_var_enc = 'Configuration';
  prf_var = 'PRF';
  prf_var_enc = 'PRF';
  ram_var = 'RAM_Taper';
  ram_var_enc = 'RAMZ20Taper';
  ram_amp_var = 'RAM';
  ram_amp_var_enc = 'RAM';
  phase_var = 'Phase_Offset_deg';
  phase_var_enc = 'PhaseZ20OffsetZ20Z28degZ29';
  wave_var_enc = 'Z23Wave';
  ttl_start_var_enc = 'TTLZ20Start';
  ttl_length_var_enc = 'TTLZ20Length';
  if ~exist('xml_file_prefix','var')
    xml_file_prefix = 'DDS';
  end
elseif xml_version == 1.2
  config_var = 'Configuration';
  config_var_enc = 'Configuration';
  prf_var = 'PRF';
  prf_var_enc = 'PRF';
  ram_var = 'RAM_Taper';
  ram_var_enc = 'RAMZ20Taper';
  ram_amp_var = 'RAM';
  ram_amp_var_enc = 'RAM';
  phase_var = 'Phase_Offset_deg';
  phase_var_enc = 'PhaseZ20OffsetZ20Z28degZ29';
  wave_var_enc = 'Z23Wave';
  ttl_start_var_enc = 'TTLZ20StartZ20Z28TTLZ30Z2DZ3ETTLZ37Z29';
  ttl_length_var_enc = 'TTLZ20LengthZ20Z28TTLZ30Z2DZ3ETTLZ37Z29';
  if ~exist('xml_file_prefix','var')
    xml_file_prefix = 'DDS';
  end
elseif xml_version == 1.3
  config_var = 'Configuration';
  config_var_enc = 'Configuration';
  prf_var = 'PRF_Hz';
  prf_var_enc = 'PRFZ20Z28HzZ29';
  ram_var = 'RAM_Taper';
  ram_var_enc = 'RAMZ20Taper';
  ram_amp_var = 'RAM';
  ram_amp_var_enc = 'RAM';
  phase_var = 'Phase_Offset_deg';
  phase_var_enc = 'PhaseZ20OffsetZ20Z28degZ29';
  wave_var_enc = 'Z23Wave';
  ttl_start_var_enc = 'TTLZ20StartZ20Z28TTLZ30Z2DZ3ETTLZ37Z29';
  ttl_length_var_enc = 'TTLZ20LengthZ20Z28TTLZ30Z2DZ3ETTLZ37Z29';
  if ~exist('xml_file_prefix','var')
    xml_file_prefix = 'DDS';
  end
elseif xml_version == 1.4
  config_var = 'DDS_Setup';
  config_var_enc = 'DDSZ20Setup';
  prf_var = 'PRF';
  prf_var_enc = 'PRF';
  ram_var = 'RAM_Taper';
  ram_var_enc = 'RAMZ20Taper';
  ram_amp_var = 'Ram_Amplitude';
  ram_amp_var_enc = 'RamZ20Amplitude';
  phase_var = 'Phase_Offset';
  phase_var_enc = 'PhaseZ20Offset';
  wave_var_enc = 'Z23Wave';
  ttl_start_var_enc = 'TTLZ20Start';
  ttl_length_var_enc = 'TTLZ20Length';
  if ~exist('xml_file_prefix','var')
    xml_file_prefix = 'radar';
  end
elseif xml_version == 1.5
  config_var = 'DDS_Setup';
  config_var_enc = 'DDSZ20Setup';
  prf_var = 'PRF';
  ram_var = 'Ram_Amplitude';
  ram_var_enc = 'RamZ20Amplitude';
  phase_var = 'Phase_Offset';
  phase_var_enc = 'PhaseZ20Offset';
  wave_var_enc = 'Z23Wave';
  ttl_start_var_enc = 'TTLZ20Start';
  ttl_length_var_enc = 'TTLZ20Length';
  if ~exist('xml_file_prefix','var')
    xml_file_prefix = 'mcords4';
  end
elseif xml_version == 1.6
  config_var = 'DDS_Setup';
  config_var_enc = 'DDSZ20Setup';
  prf_var = 'PRF';
  ram_var = 'Ram_Amplitude';
  ram_var_enc = 'RamZ20Amplitude';
  ram_amp_var = 'Ram_Amplitude';
  ram_amp_var_enc = 'RamZ20Amplitude';
  phase_var = 'Phase_Offset';
  phase_var_enc = 'PhaseZ20Offset';
  wave_var_enc = 'Z23Wave';
  ttl_start_var_enc = 'TTLZ20Start';
  ttl_length_var_enc = 'TTLZ20Length';
  NCO_freq = 'NCO_freq';
  if ~exist('xml_file_prefix','var')
    xml_file_prefix = 'mcords5';
  end
elseif xml_version == 2.0
  config_var = 'DDS_Setup';
  config_var_enc = 'DDSZ5FSetup';
  prf_var = 'PRF';
  ram_var = 'Ram_Amplitude';
  ram_var_enc = 'RamZ20Amplitude';
  ram_amp_var = 'Ram_Amplitude';
  ram_amp_var_enc = 'RamZ20Amplitude';
  phase_var = 'Phase_Offset';
  phase_var_enc = 'PhaseZ20Offset';
  wave_var_enc = 'Z23Wave';
  ttl_start_var_enc = 'TTLZ20Start';
  ttl_length_var_enc = 'TTLZ20Length';
  NCO_freq = 'NCO_freq';
  if ~exist('xml_file_prefix','var')
    xml_file_prefix = 'mcords5';
  end
else
  error('Unsupported version');
end

return;
