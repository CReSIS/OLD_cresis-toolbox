function time_after=dealwithNaN(time_before)
% This function is to fill those NaN time points. A NaN is replaced by its
% previous time points. This will generate consequtive same time points
% wihch will be dealt with time points interpolation in another function. If
% a file contains only NaN time points, it will be regarded as invalid file
% and ignored instead of data repair.This function is called by
% create_records_icards.m.
% Note: NaNs is filled and interpolated when creating records while it is
% just masked out when creating GPS files (i.e. delete those NaN
% points).---qishi
if isnan(time_before(1))% fix the first time point if it is a NaN
  time_before(1)=time_before(find(~isnan(time_before),1));
end

while any(isnan(time_before))
  idx=find(isnan(time_before),1);
  time_before(idx)=time_before(idx-1);
end

time_after=time_before;

end