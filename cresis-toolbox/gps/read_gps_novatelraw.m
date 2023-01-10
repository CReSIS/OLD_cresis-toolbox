function gps = read_gps_novatelraw(fn, param)
% gps = read_gps_novatelraw(fn, param)
%
% Parses Novatel binary commands BESTPOS,TIME
% See Novatel OEM6 or OEM7 command manuals
%
% Example:
%
%   fn = '/scratch/metadata/2022_Antarctica_BaslerMKB/20230110/GPS_Novatel_raw_aq-field22_20230110_010949.gps'
%   gps = read_gps_novatelraw(fn);
%   gps_plot(gps);
%   datestr(epoch_to_datenum(gps.gps_time(1)));
%   gps.utc_time = gps.gps_time - utc_leap_seconds(gps.gps_time(1))
% 
% Author: John Paden
%
% See also read_gps_*.m, gps_plot.m, gps_create.m

%% Read in file

% Open file
[fid,msg] = fopen(fn,'rb');
if fid < 1
  fprintf('Could not open file %s\n', fn);
  error(msg);
end

% Read entire file (may require at lot of RAM!)
A = fread(fid,inf,'uint8=>uint8');

% Close file
fclose(fid);

%% Find Frame Syncs in Novatel Data

frame_sync_byte1 = 170; % '0xAA'
frame_sync_byte2 = 68; % '0x44'
frame_sync_byte3 = 18; % '0x12'
frame_sync_byte3_or = 19; % '0x13'

frame_start_mask = A(1:end-2)==frame_sync_byte1 & A(2:end-1)==frame_sync_byte2 & (A(3:end)==frame_sync_byte3 | A(3:end)==frame_sync_byte3_or);

frame_start_idxs = find(frame_start_mask);

%% Parse each frame

% Preallocate output arrays: Assume the longest possible length for outputs
% (only BESTPOS and TIME messages present)
gps.gps_time = zeros(1,floor(length(frame_start_idxs)/2));
gps.lat = zeros(1,floor(length(frame_start_idxs)/2));
gps.lon = zeros(1,floor(length(frame_start_idxs)/2));
gps.elev = zeros(1,floor(length(frame_start_idxs)/2));
gps.roll = zeros(1,floor(length(frame_start_idxs)/2));
gps.pitch = zeros(1,floor(length(frame_start_idxs)/2));
gps.heading = zeros(1,floor(length(frame_start_idxs)/2));

next_start_idx = frame_start_idxs(1);
out_idx = 0;
have_time = false;
for frame = 1:length(frame_start_idxs)
  %fprintf('FRAME %10d of %10d\n', frame, length(frame_start_idxs));
  
  if frame_start_idxs(frame) ~= next_start_idx
    fprintf('Unexpected frame position at byte offset %d. Skipping frame.\n', frame_start_idxs(frame));
    if frame == length(frame_start_idxs)
      break;
    end
    next_start_idx = frame_start_idxs(frame+1);
    continue;
  end
  if A(frame_start_idxs(frame)+2) == 18
    message_type = 'GPS';
    header_length = double(A(frame_start_idxs(frame)+3));
    message_ID = typecast(swapbytes(A(frame_start_idxs(frame)+(4:5))),'uint16');
    gps_week = typecast(swapbytes(A(frame_start_idxs(frame)+(14:15))),'uint16');
    gps_msec = typecast(swapbytes(A(frame_start_idxs(frame)+(16:19))),'uint32');
    message_length = double(typecast(swapbytes(A(frame_start_idxs(frame)+(8:9))),'uint16'));
  else % if A(frame_start_idxs(frame)+2) == 19
    message_type = 'INS';
    header_length = 12;
    message_length = double(A(frame_start_idxs(frame)+3));
    message_ID = typecast(swapbytes(A(frame_start_idxs(frame)+(4:5))),'uint16');
    gps_week = typecast(swapbytes(A(frame_start_idxs(frame)+(6:7))),'uint16');
    gps_msec = typecast(swapbytes(A(frame_start_idxs(frame)+(8:11))),'uint32');
  end
  
  %fprintf('%s ID: %4d Length: %4d %5d Time: %4d %6.1f\n', message_type, message_ID, header_length, message_length, gps_week, gps_msec/1000);
  
  next_start_idx = frame_start_idxs(frame)+header_length+message_length+4;
  %fprintf('%g %g %g\n', next_start_idx, frame_start_idxs(frame+1),frame_start_idxs(frame+1)-next_start_idx);
  
  if message_ID==101
    % TIME MESSAGE
    %fprintf('%4d %4d %5d %6.1f\n', header_length, message_ID, gps_week, gps_msec/1000);
    year = typecast(swapbytes(A(frame_start_idxs(frame)+header_length+(28:31))),'uint32');
    month = A(frame_start_idxs(frame)+header_length+32);
    day = A(frame_start_idxs(frame)+header_length+33);
    hour = A(frame_start_idxs(frame)+header_length+34);
    minute = A(frame_start_idxs(frame)+header_length+35);
    sec = double(typecast(swapbytes(A(frame_start_idxs(frame)+header_length+(36:39))),'uint32'))/1000;
    %fprintf('%04d_%02d_%02d %02d:%02d:%04.1f\n', year, month, day, hour, minute, sec);
    have_time = true;
    
  elseif message_ID==42
    % BESTPOS MESSAGE
    %fprintf('%4d %4d %5d %6.1f\n', header_length, message_ID, gps_week, gps_msec/1000);
    if ~have_time
      fprintf('BESTPOS without preceding TIME. Skipping frame.\n', frame_start_idxs(frame));
      continue;
    end
    lat = typecast(swapbytes(A(frame_start_idxs(frame)+header_length+(8:15))),'double');
    lon = typecast(swapbytes(A(frame_start_idxs(frame)+header_length+(16:23))),'double');
    elev = typecast(swapbytes(A(frame_start_idxs(frame)+header_length+(24:31))),'double');
    %fprintf('%8.4f %8.4f %6.1f\n', lat, lon, elev);
    out_idx = out_idx + 1;
    gps.gps_time(out_idx) = datenum([double([year,month,day,hour,minute]),sec]);
    gps.lat(out_idx) = lat;
    gps.lon(out_idx) = lon;
    gps.elev(out_idx) = elev;
    have_time = false;
  end
end

%% Prepare Output

% Only keep the records that we filled:
gps.gps_time = datenum_to_epoch(gps.gps_time(1:out_idx)); % Convert from Matlab datenum to ANSI-C standard time of seconds since Jan 1, 1970
gps.gps_time = gps.gps_time + utc_leap_seconds(gps.gps_time(1)); % Convert from UTC to GPS time
gps.lat = gps.lat(1:out_idx);
gps.lon = gps.lon(1:out_idx);
gps.elev = gps.elev(1:out_idx);
gps.roll = gps.roll(1:out_idx);
gps.pitch = gps.pitch(1:out_idx);
gps.heading = gps.heading(1:out_idx);

