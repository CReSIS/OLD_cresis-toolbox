function Header = load_header_mcrds(filename)
% Header = load_header_mcrds(filename)
%
% Loads the header from an MCRDS raw file.
%
% See also: basic_load_mcrds_hdr.m, basic_load_mcrds_time.m,
%   basic_load_mcords.m, basic_load_mcords2.m, basic_load_fmcw.m,
%   basic_load_accum.m

Header = [];   

% OPEN FILE     ======================================================
fid = fopen( filename ,'r');

% READ FILE HEADER      ==============================================
[filetype num]             = fread(fid,32,'char');
Header.FileType            = char(filetype.');
[Header.Year num]          = fread(fid,1,'int16');
[Header.SubYearVer num]    = fread(fid,1,'int8');
status = fseek(fid,5,'cof');

% READ DATA HEADER      ==============================================
[datatype num]                  = fread(fid,1,'int32');
if datatype == 0; Header.IndexHeader = ftell(fid)-4;
else; fclose(fid); disp([' : ERROR : Failure to find datatype = 0 : ' filename]); end
[Header.SampleFrequency num]            = fread(fid,1,'float64');
[prfcount num]                          = fread(fid,1,'int32');
Header.PRF = 10e6/prfcount;
[Header.PreSum num]                     = fread(fid,1,'int32');
[RxAttenMatrix num]                     = fread(fid,64,'int32');
RxAttenMatrix = permute(reshape(RxAttenMatrix,[4 2 8]),[3 2 1]);
[RxBlankMatrix num]                     = fread(fid,32,'int32');
RxBlankMatrix = permute(reshape(RxBlankMatrix,[4 1 8]),[3 2 1]);
[Calibration.Enable num]                = fread(fid,1,'int8');
[Header.NumberWaveforms num]            = fread(fid,1,'int8');
status = fseek(fid,2,'cof');
[Calibration.NumberPoints num]          = fread(fid,1,'int32');
[Calibration.StartFrequency num]        = fread(fid,1,'float64');
[Calibration.StopFrequency num]         = fread(fid,1,'float64');
[Calibration.Delay num]                 = fread(fid,1,'float64');
[Calibration.Duration num]              = fread(fid,1,'float64');

% READ WAVEFORM HEADER      ==========================================
for idx = 1:Header.NumberWaveforms
    [Header.Waveform(idx).StartFrequency num]    = fread(fid,1,'float64');
    [Header.Waveform(idx).StopFrequency num]     = fread(fid,1,'float64');
    [Header.Waveform(idx).PulseDuration num]     = fread(fid,1,'float64');
    [Waveform(idx).Calibration.Frequency num]    = fread(fid,1,'float64');
    [Waveform(idx).Calibration.Delay num]        = fread(fid,1,'float64');
    [Waveform(idx).Calibration.Duration num]     = fread(fid,1,'float64');
    [Waveform(idx).BandSelect num]        = fread(fid,1,'int8');
    [Header.Waveform(idx).ZeroPiModulation num]  = fread(fid,1,'int8');
    [Header.Waveform(idx).TxMultiplexer num]     = fread(fid,1,'int8');
    [Waveform(idx).RxMode num]            = fread(fid,1,'int8');
    [Header.Waveform(idx).TxAmpEnable num]       = fread(fid,2,'int8');
    status = fseek(fid,2,'cof');
    [Header.Waveform(idx).ModCount0 num]         = fread(fid,1,'int32');
    [Header.Waveform(idx).ModCount1 num]         = fread(fid,1,'int32');
    [Header.Waveform(idx).NumberSamples num]     = fread(fid,4,'int32');
    Header.Waveform(idx).NumberSamples = ...
                Header.Waveform(idx).NumberSamples([1 1 2 2 3 3 4 4]);
    [sampledelaycount num]                       = fread(fid,4,'int32');
    Header.Waveform(idx).SampleDelay = ...
                [sampledelaycount([1 1 2 2 3 3 4 4])-10]./10e6;
    [Header.Waveform(idx).RecordEnable num]      = fread(fid,8,'int8');
    [Waveform(idx).RecordStart num]              = fread(fid,8,'int32');      
    [Waveform(idx).RecordStop num]               = fread(fid,8,'int32');
    [BlankDelayCount num]                        = fread(fid,2,'int32');
    Header.Waveform(idx).BlankDelay = BlankDelayCount/10e6;
    
    if idx == 1; Header.Waveform(idx).NumberSamples = Header.Waveform(idx).NumberSamples - ones(size(Header.Waveform(idx).NumberSamples)); end

    index = find(Header.Waveform(idx).RecordEnable==1);
    Header.Waveform(idx).NumberSamples = Header.Waveform(idx).NumberSamples(index);
    Header.Waveform(idx).SampleDelay = Header.Waveform(idx).SampleDelay(index);
    
end

% DETERMINE INDEX LOCATIONS     ==========================================
Header.IndexRecordStart = zeros(length(find(Header.Waveform(1).RecordEnable==1)),Header.NumberWaveforms);
Header.IndexRecordStop = zeros(size(Header.IndexRecordStart));
current = 0;
channel_active = find(Header.Waveform(1).RecordEnable==1).'; % if not first available recievers
for idx_channel = 1:1:length(channel_active)
    for idx_waveform = 1:1:Header.NumberWaveforms
        Header.IndexRecordStart(idx_channel,idx_waveform) = current + 1;
        % Header.IndexRecordStop(idx_channel,idx_waveform)  = Header.IndexRecordStart(idx_channel,idx_waveform) + ...
        %             Header.Waveform(idx_waveform).NumberSamples(channel_active(idx_channel)) - 1;
        Header.IndexRecordStop(idx_channel,idx_waveform)  = Header.IndexRecordStart(idx_channel,idx_waveform) + ...
                    Header.Waveform(idx_waveform).NumberSamples(idx_channel) - 1;
        current = Header.IndexRecordStop(idx_channel,idx_waveform);
    end
end

% RECORD START LOCATION OF DATA     ==================================
[datatype num]                  = fread(fid,1,'int32');
if datatype == 2; Header.IndexData = ftell(fid)-4;
else; fclose(fid); disp([' : ERROR : Failure to find datatype = 2 : ' filename]); end
status = fseek(fid,-4,'cof'); %back up so that record count is correct

% DETERMINE NUMBER OF RECORDS       ==================================
        byte_now = ftell(fid);
        fseek(fid,0,'eof');
        byte_end = ftell(fid);
        Header.NumberRecords = floor( (byte_end - byte_now) / (4+4+4+8+4+2*Header.IndexRecordStop(end,end)) );

fclose(fid);

% EDITS TO MAKE PROCESSING EASIER
% note: some fields omited from header for lack of use reasons
for idx = 1:1:Header.NumberWaveforms
    
    if Waveform(idx).BandSelect == 0
        Header.Waveform(idx).StartFrequency = Header.Waveform(idx).StartFrequency + 120e6;
        Header.Waveform(idx).StopFrequency  = Header.Waveform(idx).StopFrequency + 120e6;
    elseif Waveform(idx).BandSelect == 1
        Header.Waveform(idx).StartFrequency = Header.Waveform(idx).StartFrequency + 420e6;
        Header.Waveform(idx).StopFrequency  = Header.Waveform(idx).StopFrequency + 420e6;
    end
    
    Header.Waveform(idx).PulseDuration = Header.Waveform(idx).PulseDuration - mod(Header.Waveform(idx).PulseDuration,0.1e-6);
    
    Header.Waveform(idx).RxAttenuation = sum(RxAttenMatrix(1,:,(Waveform(idx).RxMode+1)));
    
end
    

return
    

