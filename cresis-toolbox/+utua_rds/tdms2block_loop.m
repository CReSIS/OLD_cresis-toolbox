% Michael Christoffersen
% May 2017
% Wrapper for tdms2block, converts tdms files to block files for pickgui
%
% addpath('E:\tmp\google_drive\code')
% addpath('E:\tmp\google_drive\code\tdms')

fdir = 'E:\tmp\google_drive\radar\';
files = get_filenames(fdir,'','','.tdms');
nfiles = length(files);

%% Initialize parameters output file
param_fn = fullfile(fdir,'params.csv');
fprintf('Saving %s\n', param_fn);
fid = fopen(param_fn,'w');
fprintf(fid,'file,start_time,chirp_cf,chirp_bw,chirp_len,num_traces\n');

%% Main loop
for i=1:nfiles
   fn = files{i};
   fprintf('Loading %s\n', fn);
   % Reject if not tdms file, for some reason ~endsWith(file, '.tdms') doesn't work
%    if(~endsWith(file,'.tdms'))
%        disp(endsWith(file,'.mat'))
%        continue
%    end
%    file = [fdir,'/',file];
   tdms = TDMS_getStruct(fn);
   
%    try
      % Skip to next file if less than 10 traces recorded
      if(length(tdms.meta.time.data) < 10)
         continue
      end
      block = tdms2block(tdms);
%    catch ME
%        ME.getReport
%        continue
%    end
   
   block.x(block.x == 0) = NaN;
   block.y(block.y == 0) = NaN;
   [fn_dir,fn_name] = fileparts(fn);
   out_fn = fullfile(fn_dir,[fn_name '.mat']);
   nav_fn = fullfile(fn_dir,[fn_name '.txt']);
   
   % Write out navigation
   fprintf('  Writing %s\n', nav_fn);
   dlmwrite(nav_fn,[block.lat',block.lon',block.elev_air'],'precision',10);
   
   % Write out to params file
   fprintf(fid,'%s,%s,%d,%d,%d,%d\n',fn,block.time(1),block.chirp.cf,block.chirp.bw,block.chirp.len,block.num_trace);
   
   fprintf('  Writing %s\n', out_fn);
   parsave(out_fn,block);
end

%% Cleanup

fclose(fid);

clear fdir files nfiles file tdms block outfile outnav