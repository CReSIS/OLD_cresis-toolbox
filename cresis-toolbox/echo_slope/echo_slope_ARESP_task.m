function echo_slope_ARESP_task(mdata, param)
%ECHO_SLOPE_ARESP_TASK Summary of this function goes here
%   Detailed explanation goes here

frame = mdata.Data;
[frame_rows, frame_cols] = size(frame);
avg_horiz = param.echo_slope.avg_horiz;
smoothing_start = 1 + ceil(avg_horiz / 2);
smoothing_end = frame_cols - ceil(avg_horiz / 2);
smoothed_array = zeros(frame_rows, frame_cols);

  
%% Horizontal Smoothing 

% iterate through the columns of the data and average horizontally by the
% number of pixels defined by the avg_horiz parameter
for i = 1 : frame_cols
  if i < smoothing_start
    smoothed_array(:, i) = nanmean(frame(:, 1 : (i + avg_horiz / 2)), 2);
  elseif i > smoothing_end
    smoothed_array(:, i) = nanmean(frame(:, (i - avg_horiz / 2) : frame_cols), 2);
  else
    smoothed_array(:, i) = nanmean(frame(:, (i - avg_horiz / 2) : (i + avg_horiz / 2)), 2);
  end
end

% for i = smoothing_start : smoothing_end
%     smoothed_array(:, i) = nanmean(frame(:, (i - avg_horiz / 2) : (i + avg_horiz / 2)), 2);
% end

figure(102);clf;
imagesc(lp(smoothed_array))

colormap(1-gray(256))

keyboard;

end

