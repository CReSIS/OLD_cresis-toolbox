%% Date and Time, Version
fprintf('date_time\t%s\n', datestr(now,'yyyymmdd_HHMMSS'));
A = ver('Matlab');
fprintf('version\t%s\n', A.Release);

%% FFT speed test
start_time = tic;
fprintf('fft_start\t%g\n', toc(start_time));
A = randn(10e3,5e3) + 1i*randn(10e3,5e3);
fprintf('fft_data_creation\t%g\n', toc(start_time));
for run = 1:500
  B = fft(A);
end
fprintf('fft_done\t%g\n', toc(start_time));

%% Matrix Inversion speed test
start_time = tic;
fprintf('mat_inv_start\t%g\n', toc(start_time));
A = randn(15,15,1e5) + 1i*randn(15,15,1e5);
fprintf('mat_inv_creation\t%g\n', toc(start_time));
B = zeros(size(A));
for run = 1:40
  for rline = 1:size(A,3)
    B(:,:,rline) = inv(A(:,:,rline));
  end
end
fprintf('mat_inv_done\t%g\n', toc(start_time));

%% Component wise matrix multiplies
start_time = tic;
fprintf('mat_mult_start\t%g\n', toc(start_time));
A = randn(10e3,5e3) + 1i*randn(10e3,5e3);
B = randn(10e3,5e3) + 1i*randn(10e3,5e3);
fprintf('mat_mult_creation\t%g\n', toc(start_time));
for run = 1:300
  C = A .* B;
end
fprintf('mat_mult_done\t%g\n', toc(start_time));