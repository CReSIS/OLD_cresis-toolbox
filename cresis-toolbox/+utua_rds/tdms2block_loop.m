function block = tdms2block_loop(fn,out_fn)
% tdms2block_loop(fn)
%
% Wrapper for tdms2block, converts tdms files to block files for pickgui
%
% Author: Michael Christoffersen

%% Main loop
tdms = TDMS_getStruct(fn);

if (length(tdms.meta.time.data) < 10)
  error('File has less than 10 range lines, skipping.');
end
block = utua_rds.tdms2block(tdms);

block.x(block.x == 0) = NaN;
block.y(block.y == 0) = NaN;
[fn_dir,fn_name] = fileparts(fn);

%% CReSIS Toolbox Fields

block.fn = fn;
fname = utua_rds.fname_info_utua_rds(fn);
block.datenum = fname.datenum;

block.gps_time = datenum_to_epoch(block.time);
% Assume block.time is TAI atomic time so subtract 19 leap seconds to get GPS time
block.gps_time = block.gps_time - 19;

wf = 1;
block.wfs(wf).Tpd = abs(block.chirp.len);
block.wfs(wf).fc = block.chirp.cf;
block.wfs(wf).fs = 1 / block.dt;
block.wfs(wf).dt = block.dt;
block.wfs(wf).BW = abs(block.chirp.bw);
block.wfs(wf).f0 = block.chirp.cf - block.chirp.cf*block.chirp.bw/200;
block.wfs(wf).f1 = block.chirp.cf + block.chirp.cf*block.chirp.bw/200;
block.wfs(wf).t0 = block.twtt(1);
block.wfs(wf).time = block.twtt;
block.wfs(wf).presums = block.stack;
block.wfs(wf).num_sam = size(block.ch0,1);

%% Save output
fprintf('  Writing %s\n', out_fn);
save(out_fn,'-v7.3','-struct','block');
