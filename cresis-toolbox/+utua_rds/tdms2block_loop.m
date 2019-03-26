function block = tdms2block_loop(fn,out_fn)
% tdms2block_loop(fn)
%
% Wrapper for tdms2block, converts tdms files to block files for pickgui
%
% Author: Michael Christoffersen

%% Main loop
fprintf('Loading %s\n', fn);
tdms = TDMS_getStruct(fn);

if (length(tdms.meta.time.data) < 10)
  % Need to handle too small file
  keyboard
end
block = utua_rds.tdms2block(tdms);

block.x(block.x == 0) = NaN;
block.y(block.y == 0) = NaN;
[fn_dir,fn_name] = fileparts(fn);

% Write out to params file
%fprintf(fid,'%s,%s,%d,%d,%d,%d\n',fn,block.time(1),block.chirp.cf,block.chirp.bw,block.chirp.len,block.num_trace);

fprintf('  Writing %s\n', out_fn);
utua_rds.parsave(out_fn,block);
