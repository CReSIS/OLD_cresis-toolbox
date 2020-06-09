rows = 20;
cols = 20;
flatness = 100;
along_track_weight = 100;
down_shift = rows;
matrix = magic(max([rows cols]));
matrix = matrix(1:rows, 1:cols);
matrix = matrix ./ max(matrix) .* 5;
shift_matrix = 1;
shift_elev = 1;
for i = 1:cols
  matrix(floor(abs(sin(i)/flatness)*rows+1 - shift_matrix*(i/rows)*down_shift) + down_shift, i) = 30; % Layer
end

elevation = ((down_shift - (rows:-1:1) + 1) * shift_elev + 5) ./ 1;

flattened = echo_flatten(matrix, elevation);

figure(1);
clf;
hold on;
imagesc(matrix);
colormap(1-gray);
plot(elevation);

figure(2);
clf;
hold on;
imagesc(flattened);
colormap(1-gray);
