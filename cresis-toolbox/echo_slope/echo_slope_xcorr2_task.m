function echo_slope_xcorr2_task(mdata, param )

tiles = param.echo_slope.tiles;
frame = mdata.Data;
n = param.echo_slope.n;
rows = param.echo_slope.rows;
cols = param.echo_slope.cols;


[r1, c1] = size(frame);

corr_array = zeros(r1, c1, n);
weighted_smoothing_result = zeros(r1, c1);

for i = 1:n

  C = filter2(rot90(tiles{i}.array,2), frame, 'same');
  corr_array(:,:,i) = C;

end

[~,I2] = max(corr_array,[],3);
temp = cell2mat(tiles);
theta = [temp.slope];
slope = theta(I2);

figure(1001); clf;
imagesc(slope);
colorbar
set(gcf, 'Name', 'Guassian Slope');


figure(1002); clf;
imagesc(echo_detrend(slope));
colorbar
set(gcf, 'Name', 'Detrended Guassian Slope');


smothing_size = 5;

for i = rows/2+1:r1 - rows/2
  
  fprintf('Smoothing row %d of %d at time: %s\n', i, r1-rows/2, datestr(now, 'HH:MM:SS'))
  
  for j = cols/2+1:c1-cols/2
    
    slope_tile = slope((i - smothing_size):(i + smothing_size), ((j - smothing_size):(j + smothing_size)));
    intensity_tile = frame((i - smothing_size):(i + smothing_size), ((j - smothing_size):(j + smothing_size)));
    filter_weight = intensity_tile / sum(intensity_tile(:));
    C = slope_tile.*filter_weight;
    weighted_smoothing_result(i,j,:) = sum(C(:));
  end
end

figure(1003); clf;
imagesc(echo_detrend(weighted_smoothing_result));
colorbar
set(gcf, 'Name', 'Smoothed Guassian Slope');



keyboard

end
