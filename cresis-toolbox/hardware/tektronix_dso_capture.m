function [instrument,time,data] = tektronix_dso_capture(param)
% [instrument,time,data] = tektronix_dso_capture(param)
%
% Tektronix Digital Storage Oscilloscope (DSO) waveform capture
%
% Channel names "CH1"-"CH4", "D0"-"D15"
%   Digital channels are not fully supported... Need to SELECT:BUS[01] ON?
% Trigger names can also be "AUX", "EXT"
%
% Example:
% param = [];
% param.num_ave = 4;
% param.num_pnts = 10000000;
% param.fs = 2.5e9;
% param.trigger = 'CH1';
% param.channels = {'CH1','CH2','CH3','CH4'};
% % 500 MHz, 2.5 GSPS DSO (see Utility->IO->USB Computer->USBTMC for address)
% param.visa_address = 'USB0::0x0699::0x0401::C000466::0::INSTR';
% % 1 GHz, 5 GSPS DSO (see Utility->IO->USB Computer->USBTMC for address)
% % param.visa_address = 'USB0::0x0699::0x0401::C021536::0::INSTR';
% param.time_out_sec = 60;
% [instrument,time,data] = tektronix_dso_capture(param);
%
% Authors: John Paden, Haiji Wang

num_ave = param.num_ave;
num_pnts = param.num_pnts;
fs = param.fs;
trigger = param.trigger;
channels = param.channels;
visa_address = param.visa_address;
time_out_sec = param.time_out_sec;

%% General information:
%
% See MSO4054_Programmer_Manual.pdf
%
% Does not seem to be helpful:
% out = instrhwinfo('visa')
%
% delete(dso_obj);
%
% Example of how to query DSO:
% fprintf(dso_obj, 'DATA:START?');
% start_idx = str2double(fgets(dso_obj))
%
% Examples at:
% C:\Program Files\MATLAB\R2011B\toolbox\instrument\instrument\drivers

valid_num_pnts = [1e3 1e4 1e5 1e6 1e7];
if all(num_pnts ~= valid_num_pnts)
  error('Number of points is not valid');
end

valid_num_ave = [1 2 4 8 16 32 64 128 256 512];
if all(num_ave ~= valid_num_ave)
  error('Number of averages might need to be a power of 2 (1 to 512)');
end

try
  % Create a VISA-USB object.
  % 500 MHz, 2.5 GSPS DSO (see Utility->IO->USB Computer->USBTMC for address)
  % dso_obj = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x0699::0x0401::C000466::0::INSTR', 'Tag', '');
  % 1 GHz, 5 GSPS DSO (see Utility->IO->USB Computer->USBTMC for address)
  %dso_obj = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x0699::0x0401::C021536::0::INSTR', 'Tag', '');
  
  dso_obj = instrfind('Type', 'visa-usb', 'RsrcName', param.visa_address, 'Tag', '');
  
  % Create the VISA-USB object if it does not exist otherwise use the object that was found.
  if isempty(dso_obj)
    %dso_obj = visa('NI', 'USB0::0x0699::0x0401::C000466::0::INSTR');
    %dso_obj = visa('NI', 'USB0::0x0699::0x0401::C021536::0::INSTR');
    dso_obj = visa('NI', param.visa_address);
  else
    fclose(dso_obj);
  end
  
  dso_obj = dso_obj(1);
  set(dso_obj,'InputBufferSize',num_pnts*2 + 2048);
  set(dso_obj,'Timeout',time_out_sec);
 
  fopen(dso_obj);
  
  fprintf(dso_obj, 'WFMOutpre:BYT_Nr 2');
  fprintf(dso_obj, 'DATA:ENCDG SRIBINARY');

  if num_ave == 1
    fprintf(dso_obj, ':ACQUIRE:MODE SAMPLE');
  else
    fprintf(dso_obj, ':ACQUIRE:MODE AVERAGE');
    fprintf(dso_obj, sprintf(':ACQUIRE:NUMAVG %d', num_ave));
  end
  
  fprintf(dso_obj, sprintf(':HORIZONTAL:RECORDLENGTH %d', num_pnts));
  fprintf(dso_obj, sprintf(':HORIZONTAL:SCALE %g', num_pnts/fs/10));
  
  fprintf(dso_obj, sprintf(':TRIGGER:A:EDGE:SOURCE %s', trigger));
%   fprintf(dso_obj, sprintf(':TRIGGER:A:EDGE:SLOPE %s', trigger_slope));
%   fprintf(dso_obj, sprintf(':TRIGGER:A:EDGE:LEVEL %g', trigger_level));

  for chan_idx = 1:length(channels)
    fprintf(dso_obj, 'SELECT:%s ON', channels{chan_idx});
  end
  
  fprintf(dso_obj, ':ACQUIRE:STATE OFF');
  fprintf(dso_obj, ':ACQUIRE:STOPAFTER SEQUENCE');
  fprintf(dso_obj, ':ACQUIRE:STATE ON');
  
  fprintf(dso_obj, '*OPC?');
  wait_state = fgets(dso_obj);
  if wait_state(1) ~= '1'
    error('Acquisition failed');
  end
    
  data = [];
  for chan_idx = 1:length(channels)
    channel = channels{chan_idx};
    fprintf(dso_obj,['DATA:SOURCE ' channel]);
    
    % Data Start is 1 indexed
    fprintf(dso_obj,sprintf('DATA:START %d', 1));
    fprintf(dso_obj,sprintf('DATA:STOP %d', num_pnts));

    % Clear the event status register
    fprintf(dso_obj, '*ESR?');
    fgets(dso_obj);
    
    fprintf(dso_obj, 'WAVFRM?');
    
    % Get the waveform preamble (everything up to '#')
    buffer = ''; done = false;
    while ~done
      [buffer(end+1) count] = fread(dso_obj,1,'char');
      if count == 0
        error('Early termination of waveform return');
      elseif buffer(end) == '#'
        done = true;
      end
    end
    
    fields = {};
    while ~isempty(buffer)
      [A,pos] = textscan(buffer,'%s',1,'Delimiter','#;');
      buffer = buffer(pos+1:end);
      fields{end+1} = A{1}{1};
    end
    instrument.data{chan_idx} = fields;
    
    % Get size of transfer size field
    [num_pnts_tx_size,count] = fread(dso_obj,1,'int8');
    num_pnts_tx_size = str2double(char(num_pnts_tx_size));
    
    % Get the transfer size field
    [num_bytes_tx,count] = fread(dso_obj,num_pnts_tx_size,'int8');
    num_bytes_tx = str2double(char(num_bytes_tx'));

    % Transfer data
    [data_tmp,count] = fread(dso_obj,num_bytes_tx/2,'int16');

    quantization_offset = str2double(fields{15});
    quantization_to_volts = str2double(fields{14});
    data_tmp = (data_tmp - quantization_offset) * quantization_to_volts;
    dt = str2double(fields{10});
    t0 = str2double(fields{11});
    time = t0 + dt * (0:num_pnts-1);
    
    if count == num_bytes_tx/2
      data = cat(2,data,data_tmp);
    else
      error('Number of samples read is short %d < %d', count, num_bytes_tx/2);
    end
   
    % End of line byte
    [tmp,count] = fread(dso_obj,1,'int8');
    
    % Check for errors
    fprintf(dso_obj, '*ESR?');
    instrument.esr= mat2cell(dec2bin(str2double(fgets(dso_obj)),8),1,ones(1,8));
    [instrument.PON_power_one,instrument.URQ_user_request, ...
      instrument.CME_command_error, instrument.EXE_execution_error, ...
      instrument.DDE_device_error,instrument.QYE_query_error, ...
      instrument.RQC_request_control,instrument.OPC_operation_complete] ...
      = instrument.esr{:};
    
    fprintf(dso_obj, 'EVENT?');
    instrument.event = fgets(dso_obj);
    
  end
  
%   fprintf(dso_obj,'HORIZONTAL:MAIN:SAMPLERATE?');
%   fs = str2double(fgets(dso_obj));
%   fprintf(dso_obj,'HORIZONTAL:SAMPLERATE?');
%   fs = str2double(fgets(dso_obj));
  
  fprintf(dso_obj,'*IDN?');
  instrument.idn = fgets(dso_obj);
  
catch ME
  keyboard
end

delete(dso_obj);

end
