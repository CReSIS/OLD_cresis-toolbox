function lagaan = coreg_1D(a1,a2, ident)
laggan = NaN;

% if length(a1) ~= length(a2)
%   lagaan = NaN;
%   return;
% end

debug_fig = 1;

method = {'none', 'biased', 'unbiased' , 'normalized'};

if debug_fig
  h_fig_coreg_1D = figure('Name', 'coreg_1D');
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
  legend(method, 'Location', 'east');

  subplot(121);
  plot([1:length(a2)]+lagaan, a2, 'o-');
  legend('a1','a2','a2 coreg', 'Location', 'northwest');

  %save
  [xo_table_tag, idx_xo, reuse_loc, xo_hdr] = de_ident(ident);
  sgtitle(xo_hdr, 'Interpreter', 'None');

  fig_fn = fullfile(reuse_loc, sprintf('coreg_1D_%s.png', ident));
  fprintf('Saving %s\n', fig_fn);
  ct_saveas(h_fig_coreg_1D,fig_fn);
  close(h_fig_coreg_1D)

end
