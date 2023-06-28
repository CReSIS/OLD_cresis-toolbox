% script basic_PRF_noise_analysis
%
% This script forms the outline of how to estimate the PRF of RFI from the MCoRDS
% radar data. It requires that the data be taken with no presums. We used
% it to show the noise we were seeing was 3 kHz and this turned out to be
% the ATM PRF.
%
% Instructions:
% First run basic_noise_analysis to load data.
%
% Author: John Paden

format long
if 0
  figure(1); clf;
  for rline=1:size(data,2)
    plot(lp(data(:,rline,12)))
    rline
    pause
  end
elseif 1
  figure(1); clf;
  for rline=size(data,2):-1:1
    plot(lp(data(:,rline,12)))
    rline
    pause
  end
elseif 1
  bins = [287 3054 5821 198 2965 5732 108 2875 5642 3837 1070 3926 1161]
  rlines = [4 5 6 12 13 14 20 21 22 3909 3908 3901 3900]
  clear dataf
  for bin_idx = 1:length(bins)
    dataf(:,bin_idx) = data(bins(bin_idx)+(-5:27),rlines(bin_idx),12);
  end
  imagesc(lp(abs(fft(dataf)).^2));
  plot(lp(mean(abs(fft(dataf)).^2,2)));
  
elseif 1
  bins = [287 3054 5821 198 2965 5732 108 2875 5642]
  rlines = [4 5 6 12 13 14 20 21 22]
  diff(bins)
  diff(rlines)
  % 2767 is the offset from line to line
  PRF = 9227;
  presums = 9;
  EPRF =PRF/(7+2)
  fs = 1e9/9
  noise_prf_max = fs/size(data,1)
  PRI = 1/PRF;
  EPRI = 1/EPRF
  noise_PRI_mult(4) = 1*EPRI + 2767/fs
  %noise_prf_min = 1/noise_PRI_mult(1)
  noise_PRF_factor = 1./noise_PRI_mult
  stem(noise_PRF_factor(1) * (1:18),ones(1,18))
  hold on;
  stem(noise_PRF_factor(2) * (1:18*3),ones(1,18*3),'rx')
  stem(noise_PRF_factor(3) * (1:18*5),ones(1,18*5),'c.')
  stem(noise_PRF_factor(4) * (1:18*1),ones(1,18*1),'k*')
  hold off;
elseif 1
  bins = [5652 2487 5032 1865 4412 1245 3792 625 3172]
  rlines = [3 8 11 16 19 24 27 32 35]
  diff(bins)
  diff(rlines)
  % 2547 is the offset from line to line
  PRF = 9721;
  presums = 9;
  EPRF =PRF/(7+2)
  fs = 1e9/9
  noise_prf_max = fs/size(data,1)
  PRI = 1/PRF;
  EPRI = 1/EPRF
  noise_PRI_mult(3) = 5*EPRI - 3167/fs
  %noise_prf_min = 1/noise_PRI_mult(1)
  noise_PRF_factor = 1./noise_PRI_mult
  stem(noise_PRF_factor(1) * (1:18),ones(1,18))
  hold on;
  stem(noise_PRF_factor(2) * (1:18*3),ones(1,18*3),'rx')
  stem(noise_PRF_factor(3) * (1:18*5),ones(1,18*5),'c.')
  hold off;
elseif 1
  bins = [5652 2487 5032 1865 4412 1245 3792 625 3172]
  rlines = [3 8 11 16 19 24 27 32 35]
  diff(bins)
  diff(rlines)
  % 2547 is the offset from line to line
  PRF = 9721;
  presums = 9;
  EPRF =PRF/(7+2)
  fs = 1e9/9
  noise_prf_max = fs/size(data,1)
  PRI = 1/PRF;
  EPRI = 1/EPRF
  noise_PRI_mult(2) = 3*EPRI + 2547/fs
  %noise_prf_min = 1/noise_PRI_mult(1)
  noise_PRF_factor = 1./noise_PRI_mult
  stem(noise_PRF_factor(1) * (1:18),ones(1,18))
  hold on;
  stem(noise_PRF_factor(2) * (1:18*3),ones(1,18*3),'rx')
  hold off;
elseif 1
  rlines = [1 8 9 16 17 24 25 32 33 34];
  bins = [5805 2448 5142 1782 4477 1116 3811 452 3147 5842];
  % 2695 is the offset from line to line
  PRF = 9221;
  presums = 9;
  EPRF =PRF/(7+2)
  fs = 1e9/9
  noise_prf_max = fs/size(data,1)
  PRI = 1/PRF;
  EPRI = 1/EPRF
  noise_PRI_mult(1) = EPRI + 2695/fs
  noise_prf_min = 1/noise_PRI_mult(1)
elseif 0
  % Shows outer elements are similar strength to inner elements
  colors = {'k','r','y','g','b','m','r:'};
  rline = 9;
  clear h h_label;
  figure(1); clf;
  adcs=2:16;
  for adc_idx = 1:length(adcs)
    adc = adcs(adc_idx);
    board = floor((adc-1)/4) + 1;
    rec = rline + hdrs(1).epri(1) - hdrs(board).epri(1)
    h(adc_idx) = plot(lp(data(:,rec,adc)),colors{mod(adc,length(colors))+1});
    if adc>8
      set(h(adc_idx),'LineWidth',2)
    end
    hold on;
    h_label{adc_idx} = sprintf('adc %i',adc);
  end
  legend(h,h_label)
  hold off;
end
