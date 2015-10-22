% scripts tomography_interp_poly
%
% Script for finding polynomial for tomography_interp.m
%
% Author: John Paden
%
% See also: tomography_interp_poly.m

% Try to fit to polynomial of the form:
%  [1 + X + Y + Z + X^2 + Y^2 + Z^2 + X*Y + X*Z + Y*Z + ...] * C
%  where C is a column vector of coefficients


X=0; Y=0; Z=0;
polynomial = '[1 X Y Z X^2 Y^2 X*Y X*Z Y*Z X^3 Y^3 X^2*Y Y^2*X Y^2*Z exp(X) exp(Y) X*Y^4 Y^4 Z*Y^4 Y^4  Y^5*Z Y^5]';
A = zeros(numel(value(~isnan(value))),length(eval(polynomial)));
zpos = zeros(numel(value(~isnan(value))),1);
out_idx = 0;
for row = 1:size(value,1)
  for col = 1:size(value,2)
    for third = 1:size(value,3)
      if ~isnan(value(row,col,third))
        out_idx = out_idx + 1;
        X = X_axis(row);
        Y = Y_axis(col);
        Z = Z_axis(third);
        A(out_idx,:) = [1 X Y Z X^2 Y^2 X*Y X*Z Y*Z X^3 Y^3 X^2*Y Y^2*X Y^2*Z exp(X) exp(Y) X*Y^4 Y^4 Z*Y^4 Y^4 Y^5*Z Y^5];
        zpos(out_idx) = value(row,col,third);
      end
    end
  end
end
C = pinv(A)*zpos;
round(C)

% Reconstruct the matrix
value_est = zeros(size(value));
for row = 1:size(value,1)
  for col = 1:size(value,2)
    for third = 1:size(value,3)
      if ~isnan(value(row,col,third))
        out_idx = out_idx + 1;
        X = X_axis(row);
        Y = Y_axis(col);
        Z = Z_axis(third);
        value_est(row,col,third) = [1 X Y Z X^2 Y^2 X*Y X*Z Y*Z X^3 Y^3 X^2*Y Y^2*X Y^2*Z exp(X) exp(Y) X*Y^4 Y^4 Z*Y^4 Y^4 Y^5*Z Y^5] * C;
      else
        value_est(row,col,third) = NaN;
      end
    end
  end
end

value_sub = value(:,1:60,:);
value_est_sub = value_est(:,1:60,:);

sqrt(mean(abs(value_sub(~isnan(value_sub)) - value_est_sub(~isnan(value_sub))).^2))
max(abs(value_sub(~isnan(value_sub)) - value_est_sub(~isnan(value_sub))))


figure(1); clf;
imagesc(squeeze(value(30,:,:) - value_est(30,:,:)),[-2 2]);
colorbar;

figure(2); clf;
imagesc(squeeze(value(:,30,:) - value_est(:,30,:)),[-2 2]);
colorbar;

figure(3); clf;
imagesc(squeeze(value(:,:,10) - value_est(:,:,10)),[-2 2]);
colorbar;

figure(1); clf;
imagesc(squeeze(value(2,:,:) - value_est(2,:,:)));
colorbar;

figure(2); clf;
imagesc(squeeze(value(:,1,:) - value_est(:,1,:)),[-2 2]);
colorbar;

figure(3); clf;
imagesc(squeeze(value(:,:,1) - value_est(:,:,1)),[-2 2]);
colorbar;

figure(1); clf;
imagesc(squeeze(value(end-1,:,:) - value_est(end-1,:,:)));
colorbar;

figure(2); clf;
imagesc(squeeze(value(:,end,:) - value_est(:,end,:)),[-2 2]);
colorbar;

figure(3); clf;
imagesc(squeeze(value(:,:,end) - value_est(:,:,end)),[-2 2]);
colorbar;

return;
