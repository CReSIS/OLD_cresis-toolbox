function echo_slope_task(mdata, param )

tiles = param.echo_slope.tiles;
frame = mdata.Data;
n = param.echo_slope.n;
rows = param.echo_slope.rows;
cols = param.echo_slope.cols;


[r1, c1] = size(frame);

corr_array = zeros(r1, c1, n);

for i = 1:n

  C = filter2(rot90(tiles{i}.array,2), frame, 'same');
  corr_array(:,:,i) = C;

end

corr_array = detrend(corr_array);

figure(1001); clf;
theta = linspace(param.echo_slope.min_slope,param.echo_slope.max_slope,n);
[M2,I2] = max(corr_array,[],3);
slope = theta(I2);
imagesc(slope);
colorbar





keyboard

end

