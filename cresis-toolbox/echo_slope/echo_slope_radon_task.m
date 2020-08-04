function echo_slope_radon_task(mdata, param )


frame = mdata.Data;
n = param.echo_slope.n;
rows = param.echo_slope.rows;
cols = param.echo_slope.cols;
theta = param.echo_slope.theta;

[frame_rows, frame_cols] = size(frame);

% %only look at smaller portion of the frame
% frame = frame((800):(1200), (1400):(2200));
% [frame_rows, frame_cols] = size(frame);

figure(102);
imagesc(lp(frame))

colormap(1-gray(256))

corr_array = zeros(frame_rows, frame_cols, n);
weighted_smoothing_result = zeros(frame_rows, frame_cols);

for i = rows/2+1:frame_rows - rows/2
  
  fprintf('Starting row %d of %d at time: %s\n', i, frame_rows-rows/2, datestr(now, 'HH:MM:SS'))
  
  for j = cols/2+1:frame_cols-cols/2

    tile = frame((i - rows/2):(i + rows /2), ((j - cols / 2):(j + cols /2)));
    C = sum(radon(tile.', theta).^2);
    corr_array(i,j,:) = C;
  end
end

[M2,I2] = max(corr_array,[],3);
slope = -theta(I2);

figure(1001); clf;
imagesc(slope);
colorbar
set(gcf, 'Name', 'Radon Transform Slope');


figure(1002); clf;
imagesc(echo_detrend(slope));
colorbar
set(gcf, 'Name', 'Detrended Radon Transform Slope');


for i = rows/2+1:frame_rows - rows/2
  
  fprintf('Smoothing row %d of %d at time: %s\n', i, frame_rows-rows/2, datestr(now, 'HH:MM:SS'))
  
  for j = cols/2+1:frame_cols-cols/2
    
    slope_tile = slope((i - rows/2):(i + rows /2), ((j - cols / 2):(j + cols /2)));
    intensity_tile = frame((i - rows/2):(i + rows /2), ((j - cols / 2):(j + cols /2)));
    filter_weight = intensity_tile / sum(intensity_tile(:));
    C = slope_tile.*filter_weight;
    weighted_smoothing_result(i,j,:) = sum(C(:));
  end
end

keyboard

figure(1003); clf;
imagesc(echo_detrend(weighted_smoothing_result));
colorbar
set(gcf, 'Name', 'Smoothed Radon Transform Slope');
keyboard


end

