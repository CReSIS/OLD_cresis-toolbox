rows = 20;
cols = 20;
shift_matrix = .5;
shift_elev = .5;
flatness = 10;
layer_start = 0;
down_shift = rows;

% Create default matrix
matrix = magic(max([rows cols]));
matrix = matrix(1:rows, 1:cols);
matrix = matrix ./ max(max(matrix)) .* 5;

% Create borders
matrix(1, : ) = 10;
matrix(end, : ) = 10;

% Create elevation model by which to shift the layer
layer_elev = ((down_shift - (rows:-1:1) + 1) * shift_matrix + 5) ./ 1;
elev_range = layer_elev - min(layer_elev);
elev_range = interp_finite(elev_range);

% Create layer sine wave
layer_mat = zeros(size(matrix, 1) + floor(max(elev_range)*2), size(matrix, 2));
for i = 1:cols 
  layer_mat(max(floor(abs(sin(i)/flatness)*rows) + 1, 1) + round(size(layer_mat, 1)/2), i) = 30;
end

% Skew the layer sine wave based on shift_matrix
skewed_layer_mat = layer_mat;
for i = 1:cols
  skewed_layer_mat(:, i) = interp1(1:size(skewed_layer_mat, 1), skewed_layer_mat(:, i), (1:size(matrix, 1) + floor(max(elev_range)*2)) + elev_range(i)); 
end

% remove nans
skewed_layer_mat(isnan(skewed_layer_mat)) = 0;

% Crop skewed layer matrix
[r, ~] = ind2sub([size(skewed_layer_mat, 1), size(skewed_layer_mat, 2)], find(skewed_layer_mat));
skewed_layer_mat = skewed_layer_mat(min(r):max(r), :);

% Apply layer to matrix
for c = 1:cols
   for r = 1:size(skewed_layer_mat, 1)
       row = r + layer_start;
       if skewed_layer_mat(r, c) > 0 && row <= size(matrix, 1)
           matrix(row, c) = skewed_layer_mat(r, c);
       end
   end
end

% Create elevation model by which to shift the matrix back
elevation = ((down_shift - (rows:-1:1) + 1) * shift_elev + 5) ./ 1;

flattened = echo_flatten(matrix, elevation);

% Crop layer matrix
[r, ~] = ind2sub([size(layer_mat, 1), size(layer_mat, 2)], find(layer_mat));
layer_mat = layer_mat(min(r)-5:max(r)+5, :);
figure(1);
clf;
hold on;
imagesc(layer_mat);
colormap(1-gray);
title('Layer before adding elevation shift.');

figure(2);
clf;
hold on;
imagesc(matrix);
colormap(1-gray);
plot(elevation);
legend({'Plane elevation'});
title('Layer in echogram with elevation shift.');

figure(3);
clf;
hold on;
imagesc(flattened);
colormap(1-gray);
title('Echogram shifted to account for elevation.');

