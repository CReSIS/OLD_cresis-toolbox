% scripts tomography_interp
%
% Script for finding polynomial based function which converts
% time-delay to target and direction of arrival to a y,z position.
% For use with tomography_post.m
%
% The values of interest are:
%   Cy, Cz, and scale_factors
%
% To obtain new values given height (AGL), td, and inc:
%   X = td/scale_factors(1);
%   Y = inc/scale_factors(2);
%   Z = height/scale_factors(3);
%   z_pos = [1 X Y Z X^2 Y^2 X*Y X*Z Y*Z X^3 Y^3 X^2*Y Y^2*X Y^2*Z exp(X) exp(Y) X*Y^4 Y^4 Z*Y^4 Y^4 Y^5*Z Y^5]*Cy
%   y_pos = [1 Y Z X^2 Y^2 X*Y X*Z Y*Z X^3 Y^3 X^2*Y X^2*Z Y^2*X Y^2*Z Z^2*X Y^3 Y^3*X Y^4 Y^4*Z Y^5 Y^5*Z Y^5*X]*Cz
%
% Author: John Paden
%
% See also: tomography_post.m

clear
global gRadar;
inc = (0:60)/180*pi;
fc = 195e6;
num_permittivity_pnts = 101;

% ========================

tomography_interp_tstart = tic;
physical_constants;

% Create a permitivity profile, depth_base in meters, er_base is relative
% permittivity
[depth_base,er_base] = summitPerm(fc,num_permittivity_pnts,0, ...
  fullfile(gRadar.path,'/radarSimulator/profiles/'));

% Extend for air above and down to 5000 m (deepest known ice in the world)
extra_depth = (depth_base(end) + 100 : 100 : 5400).';
depth = cat(1,-10,depth_base,extra_depth); % first element to be replaced with AGL
er = cat(1,1,er_base,er_base(end) * ones(size(extra_depth)));

z_idxs = find(depth >= 0 & depth <= 5100);

% ==============================================================
% Compute values for multiple platform heights above ground level (AGL)
% ==============================================================

% height_agls: meters above the surface
height_agls = 50:100:2000;

% Preallocate memory for outputs of loop
time_delay = zeros(length(inc), length(z_idxs), length(height_agls));
cross_track = zeros(size(time_delay));

for height_agl_idx = 1:length(height_agls)
  height_agl = height_agls(height_agl_idx);
  fprintf('AGL %d of %d (%.1f sec)\n', height_agl_idx, length(height_agls), ...
    toc(tomography_interp_tstart));
  
  % Overwrite first element of depth with current agl
  depth(1) = -height_agl;
  
  % Compute propagation time and y-position for each depth and incidence
  % angle
  for idx = 1:length(z_idxs)
    if mod(idx-1,100) == 0
      fprintf('  depth %d of %d (%.1f sec)\n', idx, length(z_idxs), ...
        toc(tomography_interp_tstart));
    end
    z_idx = z_idxs(idx);
    [time_delay(:,idx,height_agl_idx),tmp,tmp,cross_track(:,idx,height_agl_idx)] ...
      = genPropTableFromPerm(depth(1:z_idx), er(1:z_idx),fc,inc);
  end
end

% ==============================================================
% We want to be able to take time delay, incidence angle, and AGL
% and compute depth (z-position) and cross-track position (y-position)
%
% So we re-interpolate the tables so that we can plot z and y position
% versus time-delay and incidence angle and AGL.
%
% Since we already have the result vs incidence angle and AGL, we just
% have to invert depth for time delay. Since time delay is a monotonically
% increasing function wrt depth, this is easy to do with interpolation.
% ==============================================================

% time_delay/cross_track are 3-D matrixes now with the following dimensions:
%   incidence angle dependence BY depth dependence BY agl dependence
% switch to
%   depth dependence BY incidence angle dependence BY agl dependence
time_delay = permute(time_delay,[2 1 3]);
cross_track = permute(cross_track,[2 1 3]);

% For each incidence angle and agl, we invert the time_delay-depth
% relationship. First, we find the range of time_delays and create
% a uniform time delay axis to interpolate on to.
time_delay_axis = linspace(min(time_delay(:)),max(time_delay(:)),101);

y_position = zeros(length(time_delay_axis),size(time_delay,2),size(time_delay,3));
for agl_idx = 1:size(time_delay,3)
  for inc_idx = 1:size(time_delay,2)
    y_position(:,inc_idx,agl_idx) = interp1(time_delay(:,inc_idx,agl_idx), ...
      cross_track(:,inc_idx,agl_idx),time_delay_axis);
  end
end

scale_factors = [max(time_delay_axis) max(inc) max(height_agls)];
X_axis = time_delay_axis/scale_factors(1);
Y_axis = inc/scale_factors(2);
Z_axis = height_agls/scale_factors(3);

value = y_position;

tomography_interp_poly_y;
Cy = C;

z_position = zeros(length(time_delay_axis),size(time_delay,2),size(time_delay,3));
for agl_idx = 1:size(time_delay,3)
  for inc_idx = 1:size(time_delay,2)
    z_position(:,inc_idx,agl_idx) = interp1(time_delay(:,inc_idx,agl_idx), ...
      depth(z_idxs),time_delay_axis);
  end
end

value = z_position;

tomography_interp_poly_z;
Cz = C;

return;

% Example of use and a comparison (assuming you have run tomography_interp)
Cy = 1.0e+03 * [-0.001660711272909
   0.090493698716517
   0.001048136114533
   0.000217252853224
  -0.541884104307672
   4.744743290064101
   0.003971414253168
   1.309756959638642
  -0.013945276404631
   0.699992786161553
  -0.032765372052905
   0.012051703640593
  -0.001483139287919
   0.424229484875972
  -0.002626573153978
   0.699992786161575
  -0.860590076181269
  -1.524460881137507
  -0.335553751187282
   0.588158785381890
   0.903908426900442
   0.053668026221286];
Cz = 1.0e+05 * [3.508380657395748
   0.078304143910234
   3.510483542206905
  -0.011025811415163
  -0.001674155930760
   1.754974750327560
  -0.000289369503542
   0.000151185549539
   0.000429567889924
  -0.000548543415327
   0.590125520756759
   0.000283533742311
  -0.013579345401646
  -0.005777378623296
   0.002278003156437
  -3.510565644678335
   0.003606277246845
   0.065714846428169
   0.004452579897103
   0.065714846427901
  -0.007257256182412
   0.045070547933172];
scale_factors = 1.0e+03 * [0.000000095286139   0.001047197551197   1.950000000000000];
 
X = time_delay_axis(20)/scale_factors(1);
Y = inc(20)/scale_factors(2);
Z = height_agls(20)/scale_factors(3);

% y_pos = [ones(size(X)) Y Z X.^2 Y.^2 X.*Y X.*Z Y.*Z X.^3 Y.^3 X.^2.*Y X.^2.*Z Y.^2.*X Y.^2.*Z Z.^2.*X Y.^3 Y.^3.*X Y.^4 Y.^4.*Z Y.^5 Y.^5.*Z Y.^5.*X]*Cy
% z_pos = [ones(size(X)) X Y Z X.^2 Y.^2 X.*Y X.*Z Y.*Z X.^3 Y.^3 X.^2.*Y Y.^2.*X Y.^2.*Z exp(X) exp(Y) X.*Y.^4 Y.^4 Z.*Y.^4 Y.^4 Y.^5.*Z Y.^5]*Cz
y_pos = [1 Y Z X^2 Y^2 X*Y X*Z Y*Z X^3 Y^3 X^2*Y X^2*Z Y^2*X Y^2*Z Z^2*X Y^3 Y^3*X Y^4 Y^4*Z Y^5 Y^5*Z Y^5*X]*Cy
z_pos = [1 X Y Z X^2 Y^2 X*Y X*Z Y*Z X^3 Y^3 X^2*Y Y^2*X Y^2*Z exp(X) exp(Y) X*Y^4 Y^4 Z*Y^4 Y^4 Y^5*Z Y^5]*Cz

% Actual value (comparison)
y_position(20,20,20)
z_position(20,20,20)


  
  
  
