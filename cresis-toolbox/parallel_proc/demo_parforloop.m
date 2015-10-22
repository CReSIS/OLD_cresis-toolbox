% script demo_parforloop.m
%
% Example of parallel for loop
%   - fft may already be parallel?, but still shows usage

clear A
tic;
for i = 1:4
  for runs = 1:100
    A{i} = fft2(randn(1e3,1e3));
  end
end
toc;

% For the parallel for-loop use 3 cores:
matlabpool(3);

% Now rerun for loop in parallel:
clear A
tic;
parfor i = 1:4
  for runs = 1:100
    A{i} = fft2(randn(1e3,1e3));
  end
end
toc;
        
