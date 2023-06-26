function [hdr,data]=load_data_X6_v2(fn,param)
base_dir = '\\samba.cresis.ku.edu\cresis\data\HF_Sounder\2022_Greenland_Vapor\';
% load_data_X6_v2('\\samba.cresis.ku.edu\cresis\data\MCoRDS\2022_Greenland_X6\data\07282022\hell2_00010.dat',[])
%read in data
fid=fopen(fn,'r');
x =  fread(fid,[8192,1000],'int16');
header = x(1:11,:);
tmps = x(4097:4107,:);
combheader = zeros(size(header,1),size(header,2)+size(tmps,2));
combheader(:, 1:2:end) = header;
combheader(:, 2:2:end) = tmps;
header = combheader;

%waveforms
[row column] = size(header);
voltage_1=zeros(4096,column/2);
voltage_2=zeros(4096,column/2);
voltage_1 = x(12:4086,:);
voltage_2 = x(4108:8182,:);

%convert int16 to uint16
mask = header < 0;
header = header + mask*65536;
%convert two uint16s to uint32;
header = header(1:2:10,:) + 65536*header(2:2:11,:);

sync = header(1,:); %should always be 1212568132 ... ?HFRD? in ascii
epri = header(2,:); %increments one per record
time = header(4,:); %time of day in binary coded decimal
hour = bitand(bitshift(time,-20),15)*10 + bitand(bitshift(time,-16),15);
min =  bitand(bitshift(time,-12),15)*10 + bitand(bitshift(time,-8),15);
sec =  bitand(bitshift(time,-4),15)*10 + bitand(bitshift(time,-0),15);
frac = header(5,:)/1e8;
seconds = 3600*hour + 60*min + sec + frac; %seconds of day including fraction

% Reset/clear hdr struct
hdr = [];
%Load hdr struct
hdr.frame_sync = sync;
hdr.epri = epri;
hdr.seconds = seconds;
hdr.fractions = frac;
hdr.wfs_1_stop_idx = 4086;
hdr.wfs_1_start_idx = 12;
hdr.wfs_2_stop_idx = 8182;
hdr.wfs_2_start_idx = 4108;
hdr.presums = 64;

% Raw data
hdr.wfs(1).num_sam = hdr.wfs_1_stop_idx - hdr.wfs_1_start_idx+1;
hdr.wfs(2).num_sam = hdr.wfs_2_stop_idx - hdr.wfs_2_start_idx+1;
 end


