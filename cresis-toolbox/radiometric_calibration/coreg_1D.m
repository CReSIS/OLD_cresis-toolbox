function lagaan = coreg_1D(a1,a2)
laggan = NaN;

% if length(a1) ~= length(a2)
%   lagaan = NaN;
%   return;
% end

debug_fig = 1;

method = {'none', 'biased', 'unbiased' , 'normalized'};

if debug_fig
  figure;
  subplot(121);
  hold on;
  plot(a1, '.-');
  plot(a2, '.-');
  subplot(122)
  hold on;
end

try
  for m = 1:length(method)
    [r{m}, lag{m}] = xcorr(a1,a2, method{m});
    [~, max_idx] = max(r{m});
    lagaa(m) = lag{m}(max_idx);
    if debug_fig
      plot(lag{m}, r{m}, '.-');
    end
  end
end

if m==1
  laggan = lagaa;
else
  lagaan = mode(lagaa);
end

if lagaan~=0
  fprintf('coreg[none,bias,unbias,norm] [%s]=%d', num2str(lagaa), lagaan);
end

if debug_fig
  xlabel('lags');
  ylabel('xcorr');
  legend(method);

  subplot(121);
  plot([1:length(a2)]+lagaan, a2, 'o-');
  legend('a1','a2','a2 coreg');
end
