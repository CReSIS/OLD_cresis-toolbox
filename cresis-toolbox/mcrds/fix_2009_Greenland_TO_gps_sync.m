% fix_2009_Greenland_TO_gps_sync.m
%
% This function should be run to correct the records.raw.radar_time
% field in the records file. After this is done, the update_records
% function should be run to resync the GPS data.
% 
% Author: John Paden

records_fn = {};

records_fn{end+1} = '/cresis/snfs1/dataproducts/csarp_support/records/rds/2009_Greenland_TO/records_20090401_03.mat';
records_fn{end+1} = '/cresis/snfs1/dataproducts/csarp_support/records/rds/2009_Greenland_TO/records_20090402_02.mat';
records_fn{end+1} = '/cresis/snfs1/dataproducts/csarp_support/records/rds/2009_Greenland_TO/records_20090409_02.mat';
records_fn{end+1} = '/cresis/snfs1/dataproducts/csarp_support/records/rds/2009_Greenland_TO/records_20090411_01.mat';

for records_idx = 1:length(records_fn)
  records = load(records_fn{records_idx});
  
  radar_clock_delta = records.raw.radar_time - records.raw.comp_time;
  clock_correction = 0.5092 / 4.758284799999999e+03;
  radar_clock_delta = radar_clock_delta - clock_correction*(records.raw.radar_time - records.raw.radar_time(1));
  radar_clock_delta = radar_clock_delta - radar_clock_delta(1);
  
  correction = round(-radar_clock_delta*4)/4;
  
  % Make sure time is monotonically increasing (i.e. correction is monotonically increasing)
  for idx=2:length(correction)
    if correction(idx) < correction(idx-1)
      correction(idx) = correction(idx-1);
    end
  end
  
  plot(-radar_clock_delta)
  hold on
  plot(correction,'r')
  hold off
  
  records.raw.radar_time = records.raw.radar_time + correction;
  records.raw.radar_time_correction = correction;
  
  new_note = sprintf('%s correction applied %s\n', mfilename(), datestr(now));
  if ~isempty(records.notes)
    records.notes = cat(2, '\n', new_note);
  end
  records.notes = cat(2, records.notes, new_note);
  
  %keyboard % For debugging
  save(records_fn{records_idx},'-struct','records');
  
end










