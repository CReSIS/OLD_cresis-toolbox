function echo_slope_radon_task(mdata, param )


frame = mdata.Data;
n = param.echo_slope.n;
rows = param.echo_slope.rows;
cols = param.echo_slope.cols;
theta = param.echo_slope.theta;

[frame_rows, frame_cols] = size(frame);

%only look at smaller portion of the frame
frame = frame((800):(1200), (1400):(2200));
[frame_rows, frame_cols] = size(frame);

figure(102);
imagesc(lp(frame))

colormap(1-gray(256))

corr_array = zeros(frame_rows, frame_cols, n);

for i = rows/2+1:frame_rows - rows/2
  
  fprintf('Starting row %d of %d at time: %s\n', i, frame_rows-rows/2, datestr(now, 'HH:MM:SS'))
  
  for j = cols/2+1:frame_cols-cols/2

    tile = frame((i - rows/2):(i + rows /2), ((j - cols / 2):(j + cols /2)));
    C = sum(radon(tile.', theta).^2);
    corr_array(i,j,:) = C;
  end
end

%corr_array = detrend(corr_array);

figure(1001); clf;
[M2,I2] = max(corr_array,[],3);
slope = -theta(I2);
imagesc(slope);
colorbar

keyboard

end

