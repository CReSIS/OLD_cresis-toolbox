function radon_test()

clc;

%tile parameters
rows = 50;
cols = 50;
max_slope = 40;
min_slope = -40;
radon_theta = min_slope:max_slope;
n = 2;
sigma_factor = 4;

%set angle of the stub
stub_theta = 30;
stub_rows = 1000;
stub_cols = 1000;
stub_sigma_factor = 3;



frame_stub = frame_stub_generator(stub_theta, stub_rows, stub_cols, stub_sigma_factor);

% figure(100)
% imagesc(frame_stub);
% title('frame\_stub');

radon_output = sum(radon(frame_stub.', radon_theta).^2);

[M, I] = max(radon_output);



fprintf('******************************\n')
fprintf('%-20s %.3f\n','Expected Slope: ', stub_theta) 
fprintf('%-20s %.3f\n', 'Correlated Slope: ', radon_theta(I))
fprintf('******************************\n')


keyboard
 


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


frame_stub = F;
  

  

  
