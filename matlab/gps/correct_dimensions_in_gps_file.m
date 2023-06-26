% script correct_dimensions_in_gps_file


gps_fns = get_filenames(fullfile(gRadar.support_path,'gps'),'gps_20','','*.mat',struct('recursive',1));

for gps_idx = 1:length(gps_fns)
  gps_fn = gps_fns{gps_idx};
  
  fprintf('Loading %s\n', gps_fn);
  gps = load(gps_fn);
  
  fixing = false;
  if size(gps.gps_time,1) > 1
    fixing = true;
    gps.gps_time = reshape(gps.gps_time,[1 length(gps.gps_time)]);
  end
  if size(gps.gps_time,1) > 1
    fixing = true;
    gps.lat = reshape(gps.lat,[1 length(gps.lat)]);
  end
  if size(gps.gps_time,1) > 1
    fixing = true;
    gps.lon = reshape(gps.lon,[1 length(gps.lon)]);
  end
  if size(gps.gps_time,1) > 1
    fixing = true;
    gps.elev = reshape(gps.elev,[1 length(gps.elev)]);
  end
  if size(gps.gps_time,1) > 1
    fixing = true;
    gps.roll = reshape(gps.roll,[1 length(gps.roll)]);
  end
  if size(gps.gps_time,1) > 1
    fixing = true;
    gps.pitch = reshape(gps.pitch,[1 length(gps.pitch)]);
  end
  if size(gps.gps_time,1) > 1
    fixing = true;
    gps.heading = reshape(gps.heading,[1 length(gps.heading)]);
  end
  
  if fixing
    fprintf('  Fixing\n');
    save(gps_fn,'-v6','-struct','gps');
  end
  
end




