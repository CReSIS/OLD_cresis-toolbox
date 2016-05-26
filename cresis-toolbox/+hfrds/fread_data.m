function [hdr,wave] = fread_data(fid,fdat,start,num)
% [hdr,wave] = hfrds.fread_data(fid,fdat,start,num)
%
% hdr: header fields for each record
% wave: recorded waveform for each record

hdr.data_valid = zeros(num,1);
hdr.fifo = zeros(num,1);
hdr.epri = zeros(num,1);
hdr.delay = zeros(num,1);
hdr.gps = zeros(num,6);
hdr.frac = zeros(num,1);
hdr.pri = zeros(num,1);
hdr.start = zeros(num,1);
hdr.stop = zeros(num,1);
hdr.words = zeros(num,1);
hdr.bshft = zeros(num,1);
hdr.dec = zeros(num,1);
hdr.pre = zeros(num,1);
hdr.nmea_frac = zeros(num,1);
hdr.nmea = zeros(num,78);
hdr.loc = zeros(num,1);
hdr.next_loc = zeros(num,1);

for index = 1:num,
  fseek(fid,121*(start+index-1),-1);
  hdr.data_valid(index) = fread(fid,1,'uint8'); %1:1
  hdr.fifo(index) = fread(fid,1,'uint16'); %2:3
  hdr.epri(index) = fread(fid,1,'uint32'); %4:7
  hdr.delay(index) = fread(fid,1,'uint16'); %2:9
  hdr.gps(index,1:6) = char(fread(fid,6,'uint8')); %6:15
  hdr.frac(index) = fread(fid,1,'uint32'); %4:19
  hdr.pri(index) = fread(fid,1,'uint16'); %2:21
  hdr.start(index) = fread(fid,1,'uint16'); %2:23
  hdr.stop(index) = fread(fid,1,'uint16'); %2:25
  hdr.words(index) = fread(fid,1,'uint16'); %2:27
  hdr.bshft(index) = fread(fid,1,'uint8'); %1:28
  hdr.dec(index) = fread(fid,1,'uint8'); %1:29
  hdr.pre(index) = fread(fid,1,'uint16'); %2:31
  hdr.nmea_frac(index) = fread(fid,1,'uint32'); %4:35
  hdr.nmea(index,1:78) = char(fread(fid,78,'uint8')); %78:113
  hdr.loc(index) = fread(fid,1,'uint32'); %4:117
  hdr.next_loc(index) = fread(fid,1,'uint32'); %4:121
  fseek(fdat,hdr.loc(index)+128,-1);
  data = fread(fdat,hdr.words(index),'uint16');
  if index == 1, wave = zeros(hdr.words(1),num); end;
  wave(1:hdr.words(index),index) = mod(data-mod((hdr.pre(index)+1)*8192/(2^hdr.bshft(index)),65536)+32768,65536)-32768;
end;
