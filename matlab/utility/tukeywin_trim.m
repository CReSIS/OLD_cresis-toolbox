function W = tukeywin_trim(N,R)
% W = tukeywin_trim(N,R)
% 
% Wrapper function to Tukey Window (tukeywin) function that removes
% outer zeros.
%
% For example:
%  tukeywin(7,0.5)
%   ans =
%          0
%     0.7500
%     1.0000
%     1.0000
%     1.0000
%     0.7500
%          0
%  tukeywin_trim(5,0.5)
%   ans =
%     0.7500
%     1.0000
%     1.0000
%     1.0000
%     0.7500

if R == 0
  W = tukeywin(N,R);
else
  W = tukeywin(N+2,R);
  W = W(2:end-1);
end

return;
