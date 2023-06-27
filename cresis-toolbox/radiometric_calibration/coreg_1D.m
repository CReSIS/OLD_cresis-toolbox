function lagaan = coreg_1D(varargin)

% 1D = signal vector w/wo axis
% passing a windowwd
% lagaan = coreg_1D(a1,a2, ident)
% lagaan = coreg_1D(a1,x1, a2,x2, ident)

lagaan = struct();
lagaan.idxs = NaN;
lagaan.t = NaN;

lagaa = NaN;
debug_fig = 1;

switch nargin
  case 3
    % lagaan = coreg_1D(a1,a2, ident)
    basic_coreg = 1;
    a1= varargin{1};
    a2= varargin{2};
    ident= varargin{3};
  case 6
    % lagaan = coreg_1D(a1,x1, a2,x2, ident)
    basic_coreg = 0;
    a1= varargin{1};
    x1= varargin{2};
    a2= varargin{3};
    x2= varargin{4};
    interp_factor = varargin{5};
    ident= varargin{6};
end

%% basic_coreg
% if length(a1) ~= length(a2)
%   lagaan = NaN;
%   return;
% end

if basic_coreg

  % trying all
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
    lagaan.idxs = lagaa; % basic
  else
    lagaan.idxs = mode(lagaa); % skip biased
  end

  if lagaan.idxs~=0
    % would not happen if max(surface_box) is accurate in caller
    fprintf('coreg[none,bias,unbias,norm] [%s]=%d', num2str(lagaa), lagaan.idxs);
  end

  if debug_fig
    xlabel('lags');
    ylabel('xcorr');
    legend(method, 'Location', 'east');

    subplot(121);
    plot([1:length(a2)]+lagaan.idxs, a2, 'o-');
    legend('a1','a2','a2 coreg', 'Location', 'northwest');

    %save
    [xo_table_tag, idx_xo, reuse_loc, xo_hdr] = de_ident(ident);
    sgtitle(xo_hdr, 'Interpreter', 'None');

    fig_fn = fullfile(reuse_loc, sprintf('coreg_1D_%s.png', ident));
    fprintf('Saving %s\n', fig_fn);
    ct_saveas(h_fig_coreg_1D,fig_fn);
    close(h_fig_coreg_1D)
  end

else
  %% interp+coreg a2 onto a1

  interp_factor = 2*round(interp_factor/2); % even

  N_x1 = length(x1);
  N_x2 = length(x2);

  dt_x1 = median(diff(x1));
  dt_x2 = median(diff(x2));
  if dt_x2 ~= dt_x1
    keyboard;
  end

  dt_diff_tol = min([dt_x1, dt_x2]) *1e-3; % /interp_factor

  if abs(dt_x1 - dt_x2) > dt_diff_tol
      warning('Check diff(dt) before coreg');
  end

X1 = interp1(1:1:N_x1, x1, 1:1/interp_factor:N_x1);
A1 = interp1(x1, a1, X1, 'spline');

X2 = interp1(1:1:N_x2, x2, 1:1/interp_factor:N_x2);
A2 = interp1(x2, a2, X2, 'spline');


% trying all
  method = {'none', 'biased', 'unbiased' , 'normalized'};
  method = {'none','biased', 'unbiased'};

  if debug_fig
    h_fig_coreg_1D = figure('Name', 'coreg_1D');
    subplot(121);
    hold on;
    plot(A1, '.-');
    plot(A2, '.-');
    subplot(122)
    hold on;
  end

  try
    for m = 1:length(method)
      [r{m}, lag{m}] = xcorr(A1,A2, method{m});
      [~, max_idx] = max(r{m});
      lagaa(m) = lag{m}(max_idx);
      if debug_fig
        plot(lag{m}, r{m}, '.-');
      end
    end
  end

  if m==1
    lagaan.idxs = lagaa; % basic
  else
    lagaan.idxs = mode(lagaa); % skip biased
  end

  if lagaan.idxs==0
      % [~,~, D] = alignsignals(A1, A2,Method="npeak", PeakNum=1);
      % lagaan = -D;
  end

  if lagaan.idxs~=0
    % would not happen if max(surface_box) is accurate in caller
    fprintf('coreg[none,bias,unbias,norm] [%s]=%d', num2str(lagaa), lagaan.idxs);
    else
      %
  end

  lagaan.t = lagaan.idxs*dt_x2;
  
  shifted_idxs = [1:length(A2)] + lagaan.idxs;

  if debug_fig
    xlabel('lags');
    ylabel('xcorr');
    legend(method, 'Location', 'east');
    grid on; zoom on;

    subplot(121);
    plot(shifted_idxs, A2, 'o-');
    legend('A1','A2','A2 coreg', 'Location', 'northwest');
    grid on; zoom on;

    %save
    [xo_table_tag, idx_xo, reuse_loc, xo_hdr] = de_ident(ident);
    sgtitle(xo_hdr, 'Interpreter', 'None');

    fig_fn = fullfile(reuse_loc, sprintf('coreg_1D_%s.png', ident));
    set(h_fig_coreg_1D, 'Position', get(0, 'Screensize'));
    fprintf('Saving %s\n', fig_fn);
    ct_saveas(h_fig_coreg_1D,fig_fn);
    close(h_fig_coreg_1D);
  end


end