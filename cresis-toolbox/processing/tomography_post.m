% script tomography_post
%
% For taking results from blake_tomo.m (wrapper script that calls
% from tomography.m) and geocoding these.
%
% Antenna information:
% 0.857 m between RX elements,
% 3.658 m between TX elements
% and 4.365 m between the RX and TX arrays.
%

format compact; format long;
tic;
physical_constants;

% =====================================================================
% User Settings
% =====================================================================

filenames = {};

base_path = '/N/dc/projects/cresis/output/mcords/2011_Antarctica_DC8/CSARP_tomo/';
for frm = [1:10]
  filenames{end+1} = sprintf('20111014_07/Data_img_01_20111014_07_%03i',frm);
end
% base_path = '/cresis/scratch2/mdce/mcords2/2011_Greenland_P3/CSARP_music/';
% filenames{end+1} = '20110413_03/Data_img_02_20110413_03_015';
% filenames{end+1} = '20110413_03/Data_img_02_20110413_03_019';
% filenames{end+1} = '20110416_02/Data_img_02_20110416_02_019';
% filenames{end+1} = '20110413_02/Data_img_02_20110413_02_002';
% filenames{end+1} = '20110413_02/Data_img_02_20110413_02_006';
% filenames{end+1} = '20110413_02/Data_img_02_20110413_02_010';
% filenames{end+1} = '20110413_02/Data_img_02_20110413_02_014';
% filenames{end+1} = '20110413_02/Data_img_02_20110413_02_018';
% filenames{end+1} = '20110413_02/Data_img_02_20110413_02_022';
% filenames{end+1} = '20110413_02/Data_img_02_20110413_02_026';
% filenames{end+1} = '20110413_02/Data_img_02_20110413_02_030';

% fc = 195e6;
fc = 193.9e6;

% Coefficients from tomography_interp.m
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

tomography_post_loaded = zeros(length(filenames),1);
for file_idxs = 1:length(filenames)

  in_fn = fullfile(base_path,[filenames{file_idxs} '.mat']);

  yAxis_Override = linspace(-1050,1050,106);

  % Sets the debugLevel which controls how much debug information is
  % printed to the console output
  debugLevel = 1;

  % =====================================================================
  % General Setup
  % =====================================================================

  if ~tomography_post_loaded(file_idxs)
    load(in_fn);
    yAxis = yAxis_Override;
    
    time = Time;
    lat = Latitude/180*pi;
    lon = Longitude/180*pi;
    elev = Elevation;
    if 1
      % elev variable contains some sudden jumps... remove
      [Bfilt,Afilt] = butter(2,0.02);
      elev_filt = filtfilt(Bfilt,Afilt,elev);
      if debugLevel >= 4
        plot(elev);
        hold on;
        plot(elev_filt,'r');
        hold off;
        return;
      end
      elev = elev_filt;
    end

    [x,y,z] = geodetic2ecef(lat,lon,elev,WGS84.ellipsoid);
    [east,north,zone] = geodetic_to_utm(lat*180/pi,lon*180/pi);

    % Create angle of arrival axis assuming snow interface
    %   See tomography_tests.m for more explanation, but negative
    %   spatial frequencies correspond to targets on the left.
    n_snow = sqrt(1.0);
    lambda = c/(n_snow*fc);
    % dy = antenna spacing (assumed to be uniform)
    
    Nsv = size(Topography,1);
    dNsv = 1/Nsv;
    spatial_freq = dNsv * [0 : floor((Nsv-1)/2), -floor(Nsv/2) : -1];
    theta = fftshift(asin(2*spatial_freq));

    % Remove angles of arrival that are outside the visible region
    %   - The spacing of the antennas was oversampled for snow
    %     interface so some of the highest spatial frequencies
    %     output by music/fft are impossible with 160 MHz wave.
    goodInc_idxs = find(imag(theta) < 1e-5 & abs(theta) ~= pi/2);
    theta = theta(goodInc_idxs);
    body_idx = Topography(goodInc_idxs,:);

    % =====================================================================
    % First threshold: any range bin greater than 500 is NaN
    body_idx2 = body_idx;
    %body_idx2(body_idx2 > 500) = NaN;

    % =====================================================================
    % Second threshold: remove any point errors
    %   median filter is [cross-track-pixels along-track-pixels]
    %   This median filter is only used to find point errors and does not
    %   effect the data except in places where point errors are found.
    bb = medfilt2(body_idx2,[5 13]);
    cc = abs(body_idx2 - bb);
    medianThreshold = 20;
    if debugLevel >=2
      figure(1);
      imagesc(cc > medianThreshold);
      figure(2);
      imagesc(body_idx2);
      keyboard
    end
    body_idx25 = body_idx2;
    if debugLevel >= 2
      fprintf('Number of pixels thresholded: %.0f\n', sum(sum(cc(60:163,:) > medianThreshold)));
      fprintf('Number of pixels total: %.0f\n', (163-60+1)*size(cc,2));
      fprintf('Percentage: %.2f%%\n', 100*sum(sum(cc(60:163,:) > medianThreshold))/((163-60+1)*size(cc,2)));
      keyboard
    end
    body_idx25(cc > medianThreshold) = bb(cc > medianThreshold);

    
    % =====================================================================
    % Only keep values that are connected to a good start pixel
    %body_idx3 = flood(double(body_idx25),round(Ntheta/2),20);
    body_idx3 = body_idx25;

    % =====================================================================
    % Apply [cross-track-pixels along-track-pixels] median filter
    body_idx4 = medfilt2(body_idx3,[3 13]);

    % =====================================================================
    % Only keep values that are connected to a good start pixel
    %body_idx5 = flood(body_idx4,round(Ntheta/2),round(size(body_idx4,2)));
    body_idx5 = body_idx4;

    if debugLevel >= 2
      imagesc(body_idx5);
      keyboard;
    end

    % =====================================================================
    % Convert range bins to actual time delay
    dt = time(2)-time(1);
    timeDelay = time(1) + (body_idx5-1)*dt;

    % =====================================================================
    % Use polynomial to geocode results
    body_idx6_y = zeros(size(timeDelay),'single');
    body_idx6_z = zeros(size(timeDelay),'single');
    for inc_idx = 1:size(timeDelay,1)
      X = timeDelay(inc_idx,:).'/scale_factors(1);
      Y = abs(theta(inc_idx)*ones(size(X)))/scale_factors(2);
      Z = (Surface.'*c/2)/scale_factors(3);

      body_idx6_y(inc_idx,:) = sign(theta(inc_idx)) .* ...
        [ones(size(X)) Y Z X.^2 Y.^2 X.*Y X.*Z Y.*Z X.^3 Y.^3 X.^2.*Y X.^2.*Z Y.^2.*X Y.^2.*Z Z.^2.*X Y.^3 Y.^3.*X Y.^4 Y.^4.*Z Y.^5 Y.^5.*Z Y.^5.*X]*Cy;
      body_idx6_z(inc_idx,:) = ...
        [ones(size(X)) X Y Z X.^2 Y.^2 X.*Y X.*Z Y.*Z X.^3 Y.^3 X.^2.*Y Y.^2.*X Y.^2.*Z exp(X) exp(Y) X.*Y.^4 Y.^4 Z.*Y.^4 Y.^4 Y.^5.*Z Y.^5]*Cz;
    end
    
    body_idx6_z_filt = filter2(ones(5,5)/(5*5), body_idx6_z);

%     body_idx6_z_filt = -body_idx6_z_filt(53:203,5:end-4);
%     body_idx6_y = body_idx6_y(53:203,5:end-4);
% 20 to 110?
    body_idx6_z_filt = -body_idx6_z_filt(6:59,5:end-4);
    body_idx6_y = body_idx6_y(6:59,5:end-4);
%     body_idx6_z_filt = -body_idx6_z_filt(20:110,5:end-4);
%     body_idx6_y = body_idx6_y(20:110,5:end-4);
%     body_idx6_z_filt = -body_idx6_z_filt(20:53,5:end-4);
%     body_idx6_y = body_idx6_y(20:53,5:end-4);
    east = east(5:end-4);
    north = north(5:end-4);
    lat = lat(5:end-4);
    lon = lon(5:end-4);
    elev = elev(5:end-4);
    surface = Surface(5:end-4);

    % =====================================================================
    % Interpolate the cross-track MUSIC result onto uniform y-Axis... the results
    % are already on a uniform x-axis.
    %   To do this, all NaN in the current y-values must be replaced. We
    %   replace them with the average of the y-value for that spatial
    %   frequency from all the other range lines.
    % =====================================================================
    if debugLevel >= 2
      % This plot shows that the body_idx6_y y-axis does not change much
      % and is nearly the same for each spatial frequency bin suggesting
      % that the NaN substition that is used will be approximately correct.
      imagesc(body_idx6_y);
    end
    body_idx6_z_int = zeros(length(yAxis),size(body_idx6_z_filt,2));
    % For each range line
    for rline = 1:size(body_idx6_z_filt,2)
      % Make body_idx6_y monotonically increasing
      cur_val = -inf;
      for doa_bin = 1:size(body_idx6_y,1)
        if body_idx6_y(doa_bin,rline) > cur_val
          cur_val = body_idx6_y(doa_bin,rline);
        else
          body_idx6_y(doa_bin,rline) = NaN;
        end
      end
      bad_idxs = find(isnan(body_idx6_y(:,rline))).';
      % For each NaN in that range line
      for bad_idx = bad_idxs
        % Take the average of all the indices in that range bin
        tmpVals = body_idx6_y(bad_idx,:);
        body_idx6_y(bad_idx,rline) = mean(tmpVals(~isnan(tmpVals)));
      end
      % Interpolate the z estimate
      good_idx_mask = logical(zeros(size(body_idx6_y(:,rline))));
      good_idx_mask(1) = 1;
      max_cross_track_idx = 1;
      for cross_track_idx = 2:length(good_idx_mask)
        if body_idx6_y(cross_track_idx,rline) > body_idx6_y(max_cross_track_idx,rline)
          max_cross_track_idx = cross_track_idx;
          good_idx_mask(cross_track_idx) = 1;
        end
      end
      try
        body_idx6_z_int(:,rline) = interp1(body_idx6_y(good_idx_mask,rline), ...
          body_idx6_z_filt(good_idx_mask,rline),yAxis,'linear','extrap');
      catch ME
        keyboard
      end
      if 0
        plot(body_idx6_z_int(:,rline));
        keyboard
      end
    end
    if 0
      imagesc(body_idx6_z_int);
      keyboard
    end
    fprintf('  MUSIC interpolation completed (%.1f sec)\n', toc);

    tomography_post_loaded(file_idxs) = true;
  end

  xAxis = sqrt(abs(east-east(1)).^2 + abs(north-north(1)).^2);

  
  % =====================================================================
  % Velocity vector (filter to 600 m)
  [B,A] = butter(2,25/600);
  angle = atan2(diff(north),diff(east));
  angle = unwrap(angle);
  angleMean = mean(angle);
  angle = angle - angleMean;
  angleFilt = filtfilt(B,A,angle);
  angle = [angle(1) angle];
  angleFilt = [angleFilt(1) angleFilt];

  yAngle = angleFilt;

  if debugLevel >= 2
  figure(1); clf;
  plot(angle*180/pi,'r');
  grid on;
  hold on;
  plot(angleFilt*180/pi);
  hold off;
  keyboard;
  end
  
  figure; clf;
  body_idx6_z_int(body_idx6_z_int > 0) = 0;
  h = surf(xAxis/1e3,yAxis,medfilt2(double(body_idx6_z_int),[5 5],'symmetric'));
  h = surf(xAxis(1:5:end)/1e3,yAxis(1:5:end),medfilt2(double(body_idx6_z_int(1:5:end,1:5:end)),[5 5],'symmetric'));
  %h = surf(filter2(ones(5,5)/25,double(body_idx6_z_int),'valid'));
%   set(h,'EdgeAlpha',0.2);
%   view(5,60);
  zlim([-2000 0]);
  xlabel('Along-track (km)')
  ylabel('Cross-track (m)')
  zlabel('Elevation (m)');

  figure; clf;
  plot(lon*180/pi,lat*180/pi);
  xlabel('Longitude (E)');
  ylabel('Latitude (N)');

  body = body_idx6_z_int;
  save(in_fn,'lat','lon','elev','east','north','zone','surface','xAxis','yAxis','yAngle','body','angleMean','-APPEND');

  if debugLevel >= 2
    figure(1); clf;
    plot(body(:,1:100:end))
    keyboard
  end

end

return;

