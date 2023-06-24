function [hdr,data]=load_data_vapor(fn,param)
base_dir = '\\samba.cresis.ku.edu\cresis\data\HF_Sounder\2022_Greenland_Vapor\';
% date_str_dat = '07272022';
% fn_start='test98_';

%Determine number of files for the segment
%in_fns{1}=get_filenames(fullfile(base_dir,'data',date_str_dat),fn_start,'','.dat');
in_fns{1}=fn;
num_files=length(in_fns{1});



stop_idx = 8192; start_idx = 9;     %Define start and stop indices for data, header is first 8 lines
start_file=2;                       %Skip first X files

% datadir = fullfile(base_dir,date_str,'hf');
% fn = fullfile(datadir,filename);

%read in data

%x = zeros(8192,(num_files-start_file+1)*1000);      %initialize matrix to read in data
x = zeros(8192,(num_files-start_file+1)*1000);      %initialize matrix to read in data


%fid = fopen(sprintf('%s%05d.dat',fullfile(base_dir,'data',date_str_dat,fn_start),index-1),'r');
%x(1:8192,((index-start_file)*1000+1):((index-start_file+1)*1000)) =  fread(fid,[8192,1000],'int16');

fid=fopen(fn,'r');
%x(1:8192,1:1000) =  fread(fid,[8192,1000],'int16');
x =  fread(fid,[8192,1000],'int16');
fclose(fid);



% if fid < 1
%   fprintf('Could not open file %s\n', fn);
%   error(msg);
% end



data = x(start_idx:stop_idx,:);
header = x(1:start_idx-1,:);

%convert int16 to uint16
mask = header < 0;
header = header + mask*65536;
%convert two uint16s to uint32;
header = header(1:2:7,:) + 65536*header(2:2:8,:);

sync = header(1,:); %should always be 1212568132 ... ?HFRD? in ascii
epri = header(2,:); %increments one per record
time = header(3,:); %time of day in binary coded decimal
hour = bitand(bitshift(time,-20),15)*10 + bitand(bitshift(time,-16),15);
min =  bitand(bitshift(time,-12),15)*10 + bitand(bitshift(time,-8),15);
sec =  bitand(bitshift(time,-4),15)*10 + bitand(bitshift(time,-0),15);
frac = header(4,:)/100e6;
seconds = 3600*hour + 60*min + sec + frac; %seconds of day including fraction

% Reset/clear hdr struct
hdr = [];
%Load hdr struct
hdr.frame_sync = sync;
hdr.epri = epri;
hdr.seconds = seconds;
hdr.fractions = frac;
hdr.stop_idx = stop_idx;
hdr.start_idx = start_idx;
hdr.presums = 64;

% Raw data
hdr.wfs(1).num_sam = hdr.stop_idx - hdr.start_idx+1;

% All waveforms have the same start
% hdr.wfs(1).t0 = hdr.start_idx / param.fs;

 end


