function echo_slope_ARESP_task(mdata, param)
%ECHO_SLOPE_ARESP_TASK Summary of this function goes here
%   Detailed explanation goes here

frame = mdata.Data;

%constants
c  = 299792458;
c_ice = c/sqrt(3.15);
dt = mean(diff(mdata.Time)); %time of pixel
pixel_width = mdata.param_csarp.csarp.sigma_x; %width of each pixel in m
pixel_height = dt*c_ice; %height of each pixel in m

%calculate # of pixels to smooth based on input params
x_pixels = ceil(param.echo_slope.avg_x / pixel_width);
x_pixels_midpoint = ceil(x_pixels / 2);
y_pixels = ceil(param.echo_slope.avg_y / pixel_height);
y_pixels_midpoint = ceil(x_pixels / 2);

[frame_rows, frame_cols] = size(frame);

row_start = 1 + ceil(x_pixels / 2);
row_end = frame_cols - ceil(x_pixels / 2);
Px = zeros(frame_rows, frame_cols);

%image of original data
figure(101);
imagesc(lp(mdata.Data))
colormap(1-gray(256))
set(gcf, 'Name', 'Original Data');
  
%% Horizontal Smoothing 

% iterate through the columns of the data and average horizontally by the
% number of pixels defined by the x_pixels parameter
for i = 1 : frame_cols
  if i < row_start
    Px(:, i) = nanmean(frame(:, 1 : (i + x_pixels_midpoint)), 2);
  elseif i > row_end
    Px(:, i) = nanmean(frame(:, (i - x_pixels_midpoint) : frame_cols), 2);
  else
    Px(:, i) = nanmean(frame(:, (i - x_pixels_midpoint) : (i + x_pixels_midpoint)), 2);
  end
end

% for i = row_start : row_end
%     Px(:, i) = nanmean(frame(:, (i - x_pixels / 2) : (i + x_pixels / 2)), 2);
% end

figure(102);clf;
imagesc(lp(Px))
colormap(1-gray(256))
set(gcf, 'Name', 'Horizontal Smoothed Data');

%% Vertical Smoothing

Pxy1 = colfilt(Px, [y_pixels 1], 'sliding', @nanmean);

figure(103);clf;
imagesc(lp(Pxy1))
colormap(1-gray(256))
set(gcf, 'Name', 'Veritcal Smoothed Data');

Pxy2 = colfilt(Px, [(y_pixels * 4) 1], 'sliding', @nanmean);
figure(104);clf;
imagesc(lp(Pxy2))
colormap(1-gray(256))
set(gcf, 'Name', 'Veritcal Smoothed Data 2');
keyboard;

end

