% script cresis_xml_mapping
%
% Inputs:
% param.radar_name
% xml_version
%
% Outputs:

if strcmpi(param.radar_name,'mcords2')
  if xml_version == 1
    config_var = 'Configuration';
    config_var_enc = 'Configuration';
    prf_var = 'PRF_Hz';
    ram_var = 'RAM Taper';
    ram_var_enc = 'RAMZ20Taper';
    xml_file_prefix = 'DDS';
    phase_var = 'Phase_Offset_deg';
    phase_var_enc = 'PhaseZ20OffsetZ20Z28degZ29';
    wave_var_enc = 'Z23Waveforms';
    ttl_start_var_enc = 'TTLZ20StartZ20Z28TTLZ30Z2DZ3ETTLZ37Z29';
    ttl_length_var_enc = 'TTLZ20LengthZ20Z28TTLZ30Z2DZ3ETTLZ37Z29';
    TTL_prog_delay = 650;
  end
  fs = 1e9/9;
  fs_sync = 1e9/18;
elseif strcmpi(param.radar_name,'mcords3')
  if xml_version == 1
    config_var = 'Configuration';
    config_var_enc = 'Configuration';
    prf_var = 'PRF_Hz';
    ram_var = 'RAM Taper';
    ram_var_enc = 'RAMZ20Taper';
    xml_file_prefix = 'DDS';
    phase_var = 'Phase_Offset_deg';
    phase_var_enc = 'PhaseZ20OffsetZ20Z28degZ29';
    wave_var_enc = 'Z23Waveforms';
    ttl_start_var_enc = 'TTLZ20StartZ20Z28TTLZ30Z2DZ3ETTLZ37Z29';
    ttl_length_var_enc = 'TTLZ20LengthZ20Z28TTLZ30Z2DZ3ETTLZ37Z29';
    TTL_prog_delay = 650;
  elseif xml_version == 3
    config_var = 'DDS_Setup';
    config_var_enc = 'DDSZ20Setup';
    prf_var = 'PRF';
    ram_var = 'Ram_Amplitude';
    ram_var_enc = 'RamZ20Amplitude';
    xml_file_prefix = 'radar';
    phase_var = 'Phase_Offset';
    phase_var_enc = 'PhaseZ20Offset';
    wave_var_enc = 'Z23Wave';
    ttl_start_var_enc = 'TTLZ20Start';
    ttl_length_var_enc = 'TTLZ20Length';
    TTL_prog_delay = 650;
  else
    error('Unsupported version');
  end
elseif strcmpi(param.radar_name,'mcords4')
  if xml_version == 2
    config_var = 'DDS_Setup';
    config_var_enc = 'DDSZ20Setup';
    prf_var = 'PRF';
    ram_var = 'Ram_Amplitude';
    ram_var_enc = 'RamZ20Amplitude';
    xml_file_prefix = 'mcords4';
    phase_var = 'Phase_Offset';
    phase_var_enc = 'PhaseZ20Offset';
    wave_var_enc = 'Z23Wave';
    ttl_start_var_enc = 'TTLZ20Start';
    ttl_length_var_enc = 'TTLZ20Length';
    TTL_prog_delay = 650;
  else
    error('Unsupported version');
  end
elseif strcmpi(param.radar_name,'mcords5')
  if xml_version == 2
    config_var = 'DDS_Setup';
    config_var_enc = 'DDSZ5FSetup';
    prf_var = 'PRF';
    ram_var = 'Ram_Amplitude';
    ram_var_enc = 'RamZ20Amplitude';
    xml_file_prefix = 'mcords5';
    phase_var = 'Phase_Offset';
    phase_var_enc = 'PhaseZ20Offset';
    wave_var_enc = 'Z23Wave';
    ttl_start_var_enc = 'TTLZ20Start';
    ttl_length_var_enc = 'TTLZ20Length';
    NCO_freq = 'NCO_freq';
    TTL_prog_delay = 650;
  elseif xml_version == 4
    config_var = 'Configuration';
    config_var_enc = 'Configuration';
    prf_var = 'PRF';
    ram_var = 'RAM';
    ram_var_enc = 'RAM';
    xml_file_prefix = 'mcords5';
    phase_var = 'Phase_Offset_deg';
    phase_var_enc = 'PhaseZ20OffsetZ20Z28degZ29';
    wave_var_enc = 'Z23Wave';
    ttl_start_var_enc = 'TTLZ20StartZ20Z28TTLZ30Z2DZ3ETTLZ37Z29';
    ttl_length_var_enc = 'TTLZ20LengthZ20Z28TTLZ30Z2DZ3ETTLZ37Z29';
    NCO_freq = 'NCO_freq';
    TTL_prog_delay = 650;
  else
    error('Unsupported version');
  end
end

