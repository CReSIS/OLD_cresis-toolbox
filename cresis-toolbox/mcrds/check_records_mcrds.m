function records = check_records_mcrds(fn)
% records = check_records_mcrds(fn)
%
% Checks for common problems in mcords records files that may require
% recreating the records files.
%
% Author: John Paden

load(fn);

records.wfs(1)
records.wfs(1).IndexRecordStart
records.wfs(1).IndexRecordStop
for wfs = 1:length(records.wfs(1).Waveform)
  records.wfs(1).Waveform(wfs)
  records.wfs(1).Waveform(wfs).TxAmpEnable.'
  records.wfs(1).Waveform(wfs).NumberSamples.'
  records.wfs(1).Waveform(wfs).SampleDelay.'
  records.wfs(1).Waveform(wfs).RecordEnable.'
  records.wfs(1).Waveform(wfs).BlankDelay.'
end
  
for idx = 2:length(records.wfs)
  fprintf('Waveform instance %d\n', idx);
  if compare_structs(records.wfs(1),records.wfs(idx),3)
    fprintf('  Different!!!\n');
  end
end

figure(1); clf;
subplot(2,1,1);
plot(records.comp_time);
title('Computer time');
subplot(2,1,2);
plot(records.radar_time);
title('Radar time');

figure(2); clf;
plot(records.lon,records.lat);
title('Lon/Lat');

figure(3); clf;
plot(records.elev);
title('Elevation (m)');

figure(4); clf;
plot(records.roll);
title('Roll');

figure(5); clf;
plot(records.pitch);
title('Pitch');

figure(6); clf;
plot(records.heading);
title('Heading');

return;
