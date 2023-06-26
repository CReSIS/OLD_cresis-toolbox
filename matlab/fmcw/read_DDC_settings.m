function fmcw_settings = read_DDC_settings(fn_dir)
% fmcw_settings = read_DDC_settings('/home/cresis1/tmp/')

% fn_dir = '/home/cresis1/tmp/';

fns = get_filenames(fn_dir,'DDC_','','.bin');

fmcw_settings.NCO_freq = [];
fmcw_settings.nyquist_zone = [];
fmcw_settings.DC_offset = [];
fmcw_settings.input_selection = [];
fmcw_settings.DDC_or_raw_select = [];
fmcw_settings.offset = [];
fmcw_settings.DDC_filter_select = [];
fmcw_settings.computer_time = [];
for fn_idx = 3:length(fns)
  fn = fns{fn_idx};
  fid = fopen(fn,'r','ieee-be');
  fseek(fid,0,1);
  end_of_file = ftell(fid);
  fseek(fid,0,-1);
  while ftell(fid) <= end_of_file-18;
    fmcw_settings.NCO_freq(end+1) = fread(fid,1,'uint16');
    fmcw_settings.nyquist_zone(end+1) = fread(fid,1,'uint8');
    fmcw_settings.offset(end+1) = fread(fid,1,'int16');
    fmcw_settings.DDC_or_raw_select(end+1) = fread(fid,1,'char');
    fmcw_settings.input_selection(end+1) = fread(fid,1,'char');
    fmcw_settings.DC_offset(end+1) = fread(fid,1,'int16');
    fmcw_settings.DDC_filter_select(end+1) = fread(fid,1,'uint8');
    fmcw_settings.computer_time(end+1) = fread(fid,1,'uint64');
  end
  % Records:
  % 9 16-bit blocks
  % 0-1: uint16: NCO_freq
  % 2: uint8: DDC_filter_select
  % 3-4: int16: offset?
  % 5: bool: DDC_or_raw_select?
  % 6: bool: input_selection?
  % 7-8: int16: DC_offset, 
  % 9: uint8: nyquist_zone
  % 10-17: uint64: computer_time (number of seconds since Jan 1, 1904)
  % 10-13: uint32: computer_time (number of seconds since Jan 1, 1904)
  
  % fmcw_settings(1,:) unsigned 16
  %   NCO frequency?
  % fmcw_settings(2,:)
  %  0x FF 00
  %  0x FF 01
  %  0x FF 02
  %  0x FF 03
  % fmcw_settings(3,:)
  %  0x 00 4C
  %  0x 01 4C
  % fmcw_settings(4,:)
  %  0x 03 00
  %  0x 03 01
  % fmcw_settings(5,:)
  %  0x 00 20
  %  0x 01 20
  %  0x 02 20
  %  0x 03 20
  % fmcw_settings(6:7,:) == 0
  % fmcw_settings(8,:) == 26829 or 0x68CD
  % fmcw_settings(9,:)
  %  Monotonically increasing, but wrapping
  
end

fmcw_settings.computer_time = datenum(1904,1,1,0,0,fmcw_settings.computer_time);

return;
plot(fmcw_settings.NCO_freq,'x')
plot(fmcw_settings.nyquist_zone)
plot(fmcw_settings.DC_offset)
plot(fmcw_settings.input_selection)
plot(fmcw_settings.DDC_or_raw_select)
plot(fmcw_settings.offset)
plot(fmcw_settings.DDC_filter_select)
plot(fmcw_settings.computer_time)



