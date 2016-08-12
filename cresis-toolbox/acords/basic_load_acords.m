function [data, time, byte_offset] = basic_load_acords(filename,param)
% function [data, time] = basic_load_acords(filename,datatype,offset,number)
%
% This function loads in ACORDS header, transmit waveform, or radar data. 
%
% Adapted from get_data.m in the folder acrods_programs_ver2 stored with
% the ACORDS raw radar data.
%
% filename = filename of ACORDS data
% param = struct controlling loading of data
%   .datatype = 0 (header), 1 (tx waveform), 2 (data)
%   .offset = byte offset to start reading data
%   .number = number of records to load in (???)
%   .file_version = 405 (single channel, 2003_Greenland_P3) or 406
%   (multi-channel, 2004_Antarctica_P3, 2005_Greenland_TO)
%
% data = file header(s) (param.datatype = 0)
%        tx waveform (param.datatype = 1)
%        radar data (param.datatype = 2)
% time = NMEA seconds 
% byte_offset = number of byte from the beginning of the file for the
%               desired datatype.
%
% Original author: Unknown
% Modified by: Logan Smith
%

% if (nargin < 2)
%    disp('get_data: you must specify a filename and datatype');
% elseif (nargin == 2)
%    offset = 0;
%    number = inf;
% elseif (nargin == 3)
%    number = inf;
% end
if ~isfield(param,'datatype')
  error('datatype field required')
end
if ~isfield(param,'offset')
  param.offset = 0;
end
if ~isfield(param,'number')
  param.number = inf;
end
if ~isfield(param,'file_version')
  param.file_version = 406;
end
if ~isfield(param,'verbose')
  param.verbose = 0;
end
datatype = param.datatype;
offset = param.offset;
number = param.number;
file_version = param.file_version;

data = [];
time = [];

fid1=fopen (filename,'r');
if fid1==-1
   error('Could not open file %s', filename);
   return
end

% header = zeros(20,1);

% search file for number of datatype records
% turn these in a structure with named fields.  return [hdr, data]
[dt bytesreturned] = fread(fid1,1,'int32'); 
seconds = fread(fid1,1,'int32');
useconds = fread(fid1,1,'int32');
header.num_samp = fread(fid1,1,'int32');
header.presums = fread(fid1,1,'int32');
header.shifts = fread(fid1,1,'int32');
header.prf = fread(fid1,1,'double');
header.f0 = fread(fid1,1,'double');
header.f1 = fread(fid1,1,'double');
header.wf_gen_clk = fread(fid1,1,'double');
header.daq_clk = fread(fid1,1,'double');
header.tpd = fread(fid1,1,'double');
header.tx_win = fread(fid1,1,'int32');
header.samp_win_delay = fread(fid1,1,'double');
header.rx_blank_end = fread(fid1,1,'double');
header.hg_blank_end = fread(fid1,1,'double');
header.low_gain_atten = fread(fid1,1,'int32');
header.high_gain_atten = fread(fid1,1,'int32');
if file_version == 406
  header.num_elem = fread(fid1,1,'int32');
  header.elem_1 = fread(fid1,1,'int32');
  header.elem_2 = fread(fid1,1,'int32');
  header.elem_3 = fread(fid1,1,'int32');
  header.elem_4 = fread(fid1,1,'int32');
end

cont = 1;
if (datatype == 0)
   numrecords = 1;
%    print_header_info(header);
else
   numrecords = 0;
end
while cont,
   % get datatype 
   [dt bytesreturned] = fread(fid1,1,'int32');
   if bytesreturned == 0,
      break
   end
   if dt == datatype,
      numrecords = numrecords + 1;
   end
   switch dt
      case 0
         if file_version == 406
           fseek(fid1,124,'cof');
         elseif file_version == 405
           fseek(fid1,104,'cof');
         end
%         disp('Found header');
      case 1
         fseek(fid1,8,'cof');
         wflen = fread(fid1,1,'int32');
         fseek(fid1,2*wflen,'cof');
%         disp('Found Waveform');
      case 2
         if file_version == 406
           fseek(fid1,(2*header.num_samp*(header.num_elem+1))*4 + 8,'cof');
         elseif file_version == 405
           fseek(fid1,(2*header.num_samp)*2 + 8,'cof');
         end
          %disp('Found Data');
      otherwise
          warning('Bad Record Type');
   end
end
if param.verbose
disp(['Found ' num2str(numrecords) ' records of datatype ' num2str(datatype)]); 
end
if offset >= numrecords,
    %disp('Offset too large');
    fclose(fid1);
    return;
else
    if (numrecords-offset) < number
       number = numrecords-offset;
    end
end

fseek(fid1,0,'bof');

time = zeros(1,number);
switch datatype
   case 0 
%       data = zeros(20,number);
   case 1 
      data = zeros(wflen,number);
   case 2 
     if file_version == 406
       data = zeros((2*header.num_samp*(header.num_elem+1)),number);
     elseif file_version == 405
       data = zeros((2*header.num_samp),number);
     end
end

%disp('Allocated data spaces');

recordnum = 0;
offset_count = 0;
byte_offset = [];
while cont,
   % get datatype 
   [dt bytesreturned] = fread(fid1,1,'int32');
   if bytesreturned == 0,
      break
   end
   switch dt
      case 0
         if datatype == 0
            recordnum = recordnum + 1;
            if recordnum > offset
                seconds = fread(fid1,1,'int32');       
                useconds = fread(fid1,1,'int32');       
                time(1,recordnum-offset) = seconds + (1e-6 * useconds);
                data(recordnum-offset).time = seconds + (1e-6 * useconds);
                data(recordnum-offset).num_samp = fread(fid1,1,'int32');
                data(recordnum-offset).presums = fread(fid1,1,'int32');
                data(recordnum-offset).shifts = fread(fid1,1,'int32');
                data(recordnum-offset).prf = fread(fid1,1,'double');
                data(recordnum-offset).f0 = fread(fid1,1,'double');
                data(recordnum-offset).f1 = fread(fid1,1,'double');
                data(recordnum-offset).wf_gen_clk = fread(fid1,1,'double');
                data(recordnum-offset).daq_clk = fread(fid1,1,'double');
                data(recordnum-offset).tpd = fread(fid1,1,'double');
                data(recordnum-offset).tx_win = fread(fid1,1,'int32');
                data(recordnum-offset).samp_win_delay = fread(fid1,1,'double');
                data(recordnum-offset).rx_blank_end = fread(fid1,1,'double');
                data(recordnum-offset).hg_blank_end = fread(fid1,1,'double');
                data(recordnum-offset).low_gain_atten = fread(fid1,1,'int32');
                data(recordnum-offset).high_gain_atten = fread(fid1,1,'int32');
                if file_version == 406
                  data(recordnum-offset).num_elem = fread(fid1,1,'int32');
                  data(recordnum-offset).elem_1 = fread(fid1,1,'int32');
                  data(recordnum-offset).elem_2 = fread(fid1,1,'int32');
                  data(recordnum-offset).elem_3 = fread(fid1,1,'int32');
                  data(recordnum-offset).elem_4 = fread(fid1,1,'int32');
                end

                 byte_offset(1,recordnum-offset) = offset_count;
                oldheader = header;
                header = data(recordnum-offset);
%                 header.num_samp = data(1,recordnum-offset);
%                 header.num_elem = data(16,recordnum-offset);
%                 if sum(oldheader ~= header),
%                    print_header_info(header,oldheader);
%                 end 
            else
               oldheader = header;
               seconds = fread(fid1,1,'int32');    
               useconds = fread(fid1,1,'int32'); 
               header.num_samp = fread(fid1,1,'int32');
               header.presums = fread(fid1,1,'int32');
               header.shifts = fread(fid1,1,'int32');
               header.prf = fread(fid1,1,'double');
               header.f0 = fread(fid1,1,'double');
               header.f1 = fread(fid1,1,'double');
               header.wf_gen_clk = fread(fid1,1,'double');
               header.daq_clk = fread(fid1,1,'double');
               header.tpd = fread(fid1,1,'double');
               header.tx_win = fread(fid1,1,'int32');
               header.samp_win_delay = fread(fid1,1,'double');
               header.rx_blank_end = fread(fid1,1,'double');
               header.hg_blank_end = fread(fid1,1,'double');
               header.low_gain_atten = fread(fid1,1,'int32');
               header.high_gain_atten = fread(fid1,1,'int32');
               if file_version == 406
                 header.num_elem = fread(fid1,1,'int32');
                 header.elem_1 = fread(fid1,1,'int32');
                 header.elem_2 = fread(fid1,1,'int32');
                 header.elem_3 = fread(fid1,1,'int32');
                 header.elem_4 = fread(fid1,1,'int32');
               end
            end
         else
            oldheader = header;
            seconds = fread(fid1,1,'int32');
            useconds = fread(fid1,1,'int32');
            header.num_samp = fread(fid1,1,'int32');
            header.presums = fread(fid1,1,'int32');
            header.shifts = fread(fid1,1,'int32');
            header.prf = fread(fid1,1,'double');
            header.f0 = fread(fid1,1,'double');
            header.f1 = fread(fid1,1,'double');
            header.wf_gen_clk = fread(fid1,1,'double');
            header.daq_clk = fread(fid1,1,'double');
            header.tpd = fread(fid1,1,'double');
            header.tx_win = fread(fid1,1,'int32');
            header.samp_win_delay = fread(fid1,1,'double');
            header.rx_blank_end = fread(fid1,1,'double');
            header.hg_blank_end = fread(fid1,1,'double');
            header.low_gain_atten = fread(fid1,1,'int32');
            header.high_gain_atten = fread(fid1,1,'int32');
            if file_version == 406
              header.num_elem = fread(fid1,1,'int32');
              header.elem_1 = fread(fid1,1,'int32');
              header.elem_2 = fread(fid1,1,'int32');
              header.elem_3 = fread(fid1,1,'int32');
              header.elem_4 = fread(fid1,1,'int32');
            end
         end
         if file_version == 406
           offset_count = offset_count + 128;
         elseif file_version == 405
           offset_count = offset_count + 108;
         end
      case 1
         if datatype == 1
            recordnum = recordnum + 1;
            if recordnum > offset
               seconds = fread(fid1,1,'int32');       
               useconds = fread(fid1,1,'int32');       
               time(1,recordnum-offset) = seconds + (1e-6 * useconds);
               wflen = fread(fid1,1,'int32');
               data(:,recordnum-offset) = fread(fid1,wflen,'int16');           
            else
                fseek(fid1,2*wflen+12,'cof');
            end
         else
            fseek(fid1,2*wflen+12,'cof');
         end
         offset_count = offset_count + 2*wflen + 16;
      case 2
         if datatype == 2
            recordnum = recordnum + 1;
            if recordnum > offset
               seconds = fread(fid1,1,'int32');       
               useconds = fread(fid1,1,'int32'); 
               offset_count = offset_count + 12; 
               time(1,recordnum-offset) = seconds + (1e-6 * useconds);
               if file_version == 406
                 data(:,recordnum-offset) = fread(fid1,(2*header.num_samp*(header.num_elem+1)),'uint32');
                 byte_offset(1,recordnum-offset) = offset_count;
                 offset_count = offset_count + (2*header.num_samp*(header.num_elem+1))*4;
               elseif file_version == 405
                 data(:,recordnum-offset) = fread(fid1,(2*header.num_samp),'uint16');
                 byte_offset(1,recordnum-offset) = offset_count;
                 offset_count = offset_count + (2*header.num_samp)*2;
               end
            else
              if file_version == 406
                fseek(fid1,(2*header.num_samp*(header.num_elem+1))*4+8,'cof');
                offset_count = offset_count + (2*header.num_samp*(header.num_elem+1))*4 +12;
              elseif file_version == 405
                fseek(fid1,(2*header.num_samp)*2+8,'cof');
                offset_count = offset_count + (2*header.num_samp)*2 + 12;
              end
            end
         else
           if file_version == 406
             fseek(fid1,(2*header.num_samp*(header.num_elem+1))*4+8,'cof');
             offset_count = offset_count + (2*header.num_samp*(header.num_elem+1))*4 + 12;
           elseif file_version == 405
             fseek(fid1,(2*header.num_samp)*2+8,'cof');
             offset_count = offset_count + (2*header.num_samp)*2 + 12;
           end
         end
     otherwise
          warning('Bad Record Type');
   end
   if (recordnum-offset) == number
      break
   end
    
    
end



fclose(fid1);
