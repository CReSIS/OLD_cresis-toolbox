function c = redblue(m, n)
% REDBLUE Red/white/blue colormap
% 
% 1st argument (m) is number of colors you want, ranging from red to blue.
% if it's even, there will be two whites [1 1 1]
% 2nd optional argument (n) is how dark you want the end of the color range to be, e.g.,
% n=0.5 means that the red end will be [(1-n)=0.5 0 0] instead of [1 0 0]
% 
% Joe MacGregor (NASA)
% Last updated: 03/25/16

if (nargin < 1)
    m                       = size(get(gcf, 'colormap'), 1);
end

if ~mod(m, 2) % even number of colors
    m2                      = m / 2;
    switch nargin
        case {0 1}
            g               = linspace(0, 1, m2)';
            b               = ones(m2, 1);
        case 2
            b               = n + linspace(0, (1 + n), m2)';
            g               = zeros(m2, 1);
            g(b > 1)        = b(b > 1) - 1;
            b(b > 1)        = 1;
    end
    c                       = [g g b; flipud([b g g])]; % rgb
else % odd number
    m2                      = ceil(m / 2);
    switch nargin
        case {0 1}
            g               = linspace(0, 1, m2)';
            b               = ones(m2, 1);
        case 2
            b               = n + linspace(0, (1 + n), m2)';
            g               = zeros(m2, 1);
            g(b > 1)        = b(b > 1) - 1;
            b(b > 1)        = 1;
    end
    c                       = [g g b; flipud([b(1:(end - 1)) repmat(g(1:(end - 1)), 1, 2)])]; % rgb
end
