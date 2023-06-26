function hdr = load_header_acords(fn)
% [hdr] = load_header_acords(fn)
%
%  INPUTS:
%    fn = filename string
%

hdr = [];

% OPEN FILE     ======================================================
fid = fopen(fn,'r');

% READ DATA HEADER      ==============================================
[datatype num]                  = fread(fid,1,'int32');
if datatype == 0;
  hdr.Indexhdr = ftell(fid)-4;
else
  fclose(fid); disp([' : ERROR : Failure to find datatype = 0 : ' fn]);
end

headers = zeros(20, 1);
seconds = fread(fid,1,'int32');
useconds = fread(fid,1,'int32');
time = seconds + 1e-6*useconds;
headers(1:3) = fread(fid, 3, 'int32');
headers(4:9) = fread(fid, 6, 'double');
headers(10) = fread(fid,1,'int32');
headers(11:13) = fread(fid, 3, 'double');
headers(14:20) = fread(fid, 7, 'int32');
% headers(14:15) = fread(fid, 2, 'int32');
nElements = headers(16)+1;

hdr.num_sam = headers(1);
hdr.num_ave = headers(2);
hdr.num_shifts = headers(3);
hdr.PRF = headers(4);
hdr.WGC_freq = headers(7);
hdr.samp_freq = headers(8);
hdr.tx_win = headers(10); % Don't know what this is used for
hdr.num_elem = headers(16)+1;
hdr.tx_atten = 0;
hdr.rx_atten = [0, 0, 0, 0];
hdr.system_delay = 0;

% Header for the Low-Gain channel
wfs(1).pulse_width = headers(9);
wfs(1).start_freq = headers(5);
wfs(1).stop_freq = headers(6);
wfs(1).sample_delay = headers(11);
wfs(1).blank_time = headers(12);
wfs(1).tx_atten = 0;
wfs(1).rx_atten = [headers(14), headers(14), headers(14), headers(14)];
wfs(1).tx_mult = 1;
wfs(1).rx_mult = [1, 1, 1, 1];
wfs(1).num_sam = headers.num_sam;
wfs(1).slots = [headers(17), headers(18), headers(19), headers(20)];
wfs(1).zero_pi_mod = 0;
wfs(1).use_ref = 0;

% Header for the High-Gain channel
wfs(2).pulse_width = headers(9);
wfs(2).start_freq = headers(5);
wfs(2).stop_freq = headers(6);
wfs(2).sample_delay = headers(11);
wfs(2).blank_time = headers(13);
wfs(2).tx_atten = 0;
wfs(2).rx_atten = [headers(15), headers(15), headers(15), headers(15)];
wfs(2).tx_mult = 1;
wfs(2).rx_mult = [1, 1, 1, 1];
wfs(2).num_sam = headers.num_sam;
wfs(2).slots = [headers(17), headers(18), headers(19), headers(20)];
wfs(2).zero_pi_mod = 0;
wfs(2).use_ref = 0;

hdr.num_wfs = 2;
hdr.dsp_mode = 0;

[hdr.fs num]            = fread(fid,1,'float64');
[prfcount num]                          = fread(fid,1,'int32');
hdr.PRF = 10e6/prfcount;
[hdr.PreSum num]                     = fread(fid,1,'int32');
[RxAttenMatrix num]                     = fread(fid,64,'int32');
RxAttenMatrix = permute(reshape(RxAttenMatrix,[4 2 8]),[3 2 1]);
[RxBlankMatrix num]                     = fread(fid,32,'int32');
RxBlankMatrix = permute(reshape(RxBlankMatrix,[4 1 8]),[3 2 1]);
[Calibration.Enable num]                = fread(fid,1,'int8');
[hdr.NumberWaveforms num]            = fread(fid,1,'int8');
status = fseek(fid,2,'cof');
[Calibration.NumberPoints num]          = fread(fid,1,'int32');
[Calibration.StartFrequency num]        = fread(fid,1,'float64');
[Calibration.StopFrequency num]         = fread(fid,1,'float64');
[Calibration.Delay num]                 = fread(fid,1,'float64');
[Calibration.Duration num]              = fread(fid,1,'float64');

% READ WAVEFORM HEADER      ==========================================
for idx = 1:hdr.NumberWaveforms
  [hdr.wfs(idx).StartFrequency num]    = fread(fid,1,'float64');
  [hdr.wfs(idx).StopFrequency num]     = fread(fid,1,'float64');
  [hdr.wfs(idx).PulseDuration num]     = fread(fid,1,'float64');
  [wfs(idx).Calibration.Frequency num]    = fread(fid,1,'float64');
  [wfs(idx).Calibration.Delay num]        = fread(fid,1,'float64');
  [wfs(idx).Calibration.Duration num]     = fread(fid,1,'float64');
  [wfs(idx).BandSelect num]        = fread(fid,1,'int8');
  [hdr.wfs(idx).ZeroPiModulation num]  = fread(fid,1,'int8');
  [hdr.wfs(idx).TxMultiplexer num]     = fread(fid,1,'int8');
  [wfs(idx).RxMode num]            = fread(fid,1,'int8');
  [hdr.wfs(idx).TxAmpEnable num]       = fread(fid,2,'int8');
  status = fseek(fid,2,'cof');
  [hdr.wfs(idx).ModCount0 num]         = fread(fid,1,'int32');
  [hdr.wfs(idx).ModCount1 num]         = fread(fid,1,'int32');
  [hdr.wfs(idx).NumberSamples num]     = fread(fid,4,'int32');
  hdr.wfs(idx).NumberSamples = ...
    hdr.wfs(idx).NumberSamples([1 1 2 2 3 3 4 4]);
  [sampledelaycount num]                       = fread(fid,4,'int32');
  hdr.wfs(idx).SampleDelay = ...
    [sampledelaycount([1 1 2 2 3 3 4 4])-10]./10e6;
  [hdr.wfs(idx).RecordEnable num]      = fread(fid,8,'int8');
  [wfs(idx).RecordStart num]              = fread(fid,8,'int32');
  [wfs(idx).RecordStop num]               = fread(fid,8,'int32');
  [BlankDelayCount num]                        = fread(fid,2,'int32');
  hdr.wfs(idx).BlankDelay = BlankDelayCount/10e6;
  
  if idx == 1; hdr.wfs(idx).NumberSamples = hdr.wfs(idx).NumberSamples - ones(size(hdr.wfs(idx).NumberSamples)); end
  
  index = find(hdr.wfs(idx).RecordEnable==1);
  hdr.wfs(idx).NumberSamples = hdr.wfs(idx).NumberSamples(index);
  hdr.wfs(idx).SampleDelay = hdr.wfs(idx).SampleDelay(index);
  
end

% DETERMINE INDEX LOCATIONS     ==========================================
hdr.IndexRecordStart = zeros(length(find(hdr.wfs(1).RecordEnable==1)),hdr.Numberwfss);
hdr.IndexRecordStop = zeros(size(hdr.IndexRecordStart));
current = 0;
channel_active = find(hdr.wfs(1).RecordEnable==1).'; % if not first available recievers
for idx_channel = 1:1:length(channel_active)
  for idx_waveform = 1:1:hdr.Numberwfss
    hdr.IndexRecordStart(idx_channel,idx_waveform) = current + 1;
    % hdr.IndexRecordStop(idx_channel,idx_waveform)  = hdr.IndexRecordStart(idx_channel,idx_waveform) + ...
    %             hdr.wfs(idx_waveform).NumberSamples(channel_active(idx_channel)) - 1;
    hdr.IndexRecordStop(idx_channel,idx_waveform)  = hdr.IndexRecordStart(idx_channel,idx_waveform) + ...
      hdr.wfs(idx_waveform).NumberSamples(idx_channel) - 1;
    current = hdr.IndexRecordStop(idx_channel,idx_waveform);
  end
end

% RECORD START LOCATION OF DATA     ==================================
[datatype num]                  = fread(fid,1,'int32');
if datatype == 2; hdr.IndexData = ftell(fid)-4;
else fclose(fid); disp([' : ERROR : Failure to find datatype = 2 : ' fn]); end
status = fseek(fid,-4,'cof'); %back up so that record count is correct

% DETERMINE NUMBER OF RECORDS       ==================================
byte_now = ftell(fid);
fseek(fid,0,'eof');
byte_end = ftell(fid);
hdr.NumberRecords = floor( (byte_end - byte_now) / (4+4+4+8+4+2*hdr.IndexRecordStop(end,end)) );

fclose(fid);

% EDITS TO MAKE PROCESSING EASIER
% note: some fields omited from header for lack of use reasons
for idx = 1:1:hdr.Numberwfss
  
  if wfs(idx).BandSelect == 0
    hdr.wfs(idx).StartFrequency = hdr.wfs(idx).StartFrequency + 120e6;
    hdr.wfs(idx).StopFrequency  = hdr.wfs(idx).StopFrequency + 120e6;
  elseif wfs(idx).BandSelect == 1
    hdr.wfs(idx).StartFrequency = hdr.wfs(idx).StartFrequency + 420e6;
    hdr.wfs(idx).StopFrequency  = hdr.wfs(idx).StopFrequency + 420e6;
  end
  
  hdr.wfs(idx).PulseDuration = hdr.wfs(idx).PulseDuration - mod(hdr.wfs(idx).PulseDuration,0.1e-6);
  
  hdr.wfs(idx).RxAttenuation = sum(RxAttenMatrix(1,:,(wfs(idx).RxMode+1)));
  
end


return


