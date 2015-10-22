function md5_copy(src,dest,qc_level)

if ~exist('qc_level','var')
  qc_level = -1;
end

while src(end) == filesep
  src = src(1:end-1);
end
src(end+1) = filesep;
fns = get_filenames(src,'','','',struct('recursive',1));

if 0
  nas = ftp('172.18.1.15','cresis1','Goosah^30');
end
bytes = 0;
old_dest_fn_dir = NaN;
for fn_idx = 1:length(fns)
  src_fn = fns{fn_idx};
  if exist(src_fn,'dir')
    continue;
  end
  src_fn_extra = src_fn(length(src)+1:end);
  fprintf('  %s %d of %d (%s)', src_fn_extra, fn_idx, length(fns), datestr(now,'HH:MM:SS'));
  dest_fn = fullfile(dest,src_fn_extra);
  dest_fn_dir = fileparts(dest_fn);
  if ~strcmpi(old_dest_fn_dir,dest_fn_dir)
    if ~exist(dest_fn_dir,'dir')
      mkdir(dest_fn_dir);
    end
  end
  
  if ~exist(dest_fn)
    done = false;
  else
    if qc_level == 1
      cmd = sprintf('md5sum %s | awk {''print $1''}',src_fn);
      [status1,result1] = system(cmd);
      cmd = sprintf('md5sum %s | awk {''print $1''}',dest_fn);
      [status2,result2] = system(cmd);
      if ~strcmp(result1,result2)
        done = false;
      else
        done = true;
      end
    else
      finfo1 = dir(src_fn);
      finfo2 = dir(dest_fn);
      if finfo1.datenum <= finfo2.datenum && finfo1.bytes == finfo2.bytes
        done = true;
      else
        done = false;
      end
    end
  end
  
  if ~done
    finfo1 = dir(src_fn);
    bytes = bytes + finfo1.bytes;
    
    while ~done
      if 1
        copyfile(src_fn,dest_fn,'f');
      else
        ftp_dest_fn_dir = ['/mnt/array1' dest_fn_dir(5:end)];
        cd(nas,ftp_dest_fn_dir);
        mput(nas,src_fn);
      end
      
      if qc_level >= 0
        cmd = sprintf('md5sum %s | awk {''print $1''}',src_fn);
        [status1,result1] = system(cmd);
        cmd = sprintf('md5sum %s | awk {''print $1''}',dest_fn);
        [status2,result2] = system(cmd);
        if ~strcmp(result1,result2)
          done = true;
        end
      else
        finfo1 = dir(src_fn);
        finfo2 = dir(dest_fn);
        done = true;
        if finfo1.datenum <= finfo2.datenum && finfo1.bytes == finfo2.bytes
          done = true;
        end
      end
      if ~done
        fprintf(' failed');
      end
    end
    fprintf(' copied (%d / %d MB total)\n', round(finfo1.bytes/2^20), ...
      round(bytes/2^20));
  else
    fprintf(' skipping\n');
  end
end

end

