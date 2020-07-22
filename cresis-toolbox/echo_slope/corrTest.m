function corrTest()

clc;

%tile parameters
rows = 100;
cols = 100;
max_slope = 10;
min_slope = -10;
n = 2;
sigma_factor = 4;

%set angle of the stub
stub_theta = 6;
stub_rows = 1000;
stub_cols = 1000;
stub_sigma_factor = 3;
norm_tiles = true;



tiles = makeTiles(rows, cols, max_slope, min_slope, n, sigma_factor, norm_tiles);


frame_stub = frame_stub_generator(stub_theta, stub_rows, stub_cols, stub_sigma_factor);

figure(100)
imagesc(frame_stub);
title('frame\_stub');



% corr_array = corr_check(tiles, frame_stub, n);
% 
% if 1
%   [M, I] = max(corr_array(round((stub_rows)/2), round((stub_cols)/2), :));
%   squeeze(corr_array(round((stub_rows)/2), round((stub_cols)/2), :))
% elseif 0
%   [M, I] = max(corr_array(round((rows+stub_rows)/2), round((cols+stub_cols)/2), :));
%   squeeze(corr_array(round((rows+stub_rows)/2), round((cols+stub_cols)/2), :))
% else
%   [M,I] = max(max(max(corr_array,[],1),[],2),[],3);
% end
% 
% figure(1001); clf;
% theta = linspace(min_slope,max_slope,n);
% [M2,I2] = max(corr_array,[],3);
% slope = theta(I2);
% imagesc(slope);
% colorbar


% figure(13)
% imagesc(corr_array(:,:,I));

% fprintf('******************************\n')
% fprintf('%-20s %.3f\n','Expected Slope: ', stub_theta) 
% fprintf('%-20s %.3f\n', 'Correlated Slope: ', tiles{I}.slope)
% fprintf('******************************\n')

% for i = 1:n
%   SUM(i) = sum(tiles{i}.array(:));
% end
% 
% figure(12)
% plot(SUM);

keyboard
 
function tiles = makeTiles(rows, cols, max_slope, min_slope, n, sigma_factor, norm)
  
%sigma X and Y based on the rows and columns
sigmaX = 1/(15 * sigma_factor);
sigmaY = 1.5/(sigma_factor);

x = -cols/2:cols/2;
y = -rows/2:rows/2;

%array of theta angles 
theta = linspace(min_slope, max_slope, n);

%cell array of tiles 
tiles = {n};

%populate each tile with slope data 
for i = 1:n
  [X,Y] = meshgrid(x,y);
  a = 12 * (pi) * sqrt(sigmaX^2 * sigmaY^2);
  B = (((X.*cosd(-theta(i)) - Y*sind(-theta(i))).^2) ./ 2*sigmaX^2);
  C = (((X.*sind(-theta(i)) - Y*cosd(-theta(i))).^2) ./ 2*sigmaY^2);
  F = a*exp(-(B+C));
  
  %normalize so sum of F = 1
  if norm
    F = F / sum(F(:));
  end
  
  tile.array = F;
  tile.slope = theta(i);

  tiles{i} = tile;

  figure(i);
  imagesc(tiles{i}.array);
  image_label = sprintf('Tile # %d, theta = %.3f', i, theta(i));
  title(image_label);
end

function frame_stub = frame_stub_generator(theta, rows, cols, sigma_factor)

sigmaX = 1/(15 * sigma_factor);
sigmaY = 1/(sigma_factor);

x = -cols/2:cols/2;
y = -rows/2:rows/2;

  
[X,Y] = meshgrid(x,y);
a = 12 * (pi) * sqrt(sigmaX^2 * sigmaY^2);
B = (((X.*cosd(theta) - Y*sind(theta)).^2) ./ 2*sigmaX^2);
C = (((X.*sind(theta) - Y*cosd(theta)).^2) ./ 2*sigmaY^2);
F = a*exp(-(B+C));
F = F / sum(F(:));

frame_stub = F;
  
function corr_array = corr_check(tiles, frameStub, n)
  
[r1, c1] = size(tiles{1}.array);
[r2, c2] = size(frameStub);
 
% corr_rows = r1 + r2 - 1;
% corr_cols = c1 + c2 - 1;
 
corr_rows = r2;
corr_cols = c2;
  
corr_array = zeros(corr_rows, corr_cols, n);

for i = 1:n

%   C = xcorr2(frameStub, tiles{i}.array);
  C = filter2(flipud(fliplr(tiles{i}.array)), frameStub, 'same');
  corr_array(:,:,i) = C;

end

  

  
