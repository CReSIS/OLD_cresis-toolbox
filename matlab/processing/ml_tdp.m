function u = ml_tdp(tx,fc,time,xx,href,hz,st,ux,alen,bins,start_range_bin,end_range_bin,er)
% function [u,tx,tau,ix0,ix1] = ml_tdp(tx,fc,time,xx,hz,st,ux,alen,bins,er)
% [u,tx,tau,ix0,ix1] = ml_tdp(tx,fc,time,xx,hz,st,ux,alen,bins,er)
%
% tx = Nt_in by Nx_in time domain data
% fc = center frequency (Hz)
% time = time axis for tx, 1 by Nt_in (seconds)
% xx = 1 by Nx_in along track vector of measurements (meters)
% hz = antenna height 1 by Nx_in (meters relative to WGS-84 ellipsoid)
% st = surface profile 1 by Nx_in (two way travel time to surface, seconds)
% ux = output image x-positions, 1 by Nx (meters)
% alen = SAR aperture length (meters)
% bins = Number of subapertures, nonnegative integer
% er = dielectric of medium to use (typical is 3.15)
%
% u = Nx by Nt by bins output image data matrix

physical_constants;

dt = time(2)-time(1);

% Compute the surface height relative to WGS-84 ellipsoid
% sz = hz - st * c/2;
sz = href - st * c/2;

%calculate local surface angle
sm = conv(sz,[1 0 -1],'same')./conv(xx,[1 0 -1],'same'); sm(1)=sm(2); sm(end)=sm(end-1);
sm = atan(sm);
sx = xx; hx = xx;

%downconvert tx
%tx = tx.*(exp(-1i*2*pi*fc*dt*(1:n(1))')*ones(1,n(2)));

%allocate output matrix
nuz = size(tx,1);
nux = numel(ux);
u = zeros(nuz,nux,bins);

su = interp1(sx,sz,ux); %local surface height at the image location (used to determine if the pixel is above or below the surface)
% hu = interp1(hx,hz,ux); %local height
hu = interp1(hx,href,ux,'linear','extrap'); %local height, with href, both hu and uz is the same for all channels

%loop over image x position
for u_idx = 1:nux %346:593 
% for u_idx = 346:593 

%     su = interp1(sx,sz,ux(u_idx)); %local surface height at the image location (used to determine if the pixel is above or below the surface)
%     hu = interp1(hx,hz,ux(u_idx)); %local height
%     agl_eff = hu-su;
%     agl_eff = hu(u_idx)-su(u_idx);

    %determine local range of x indexes for h and s
    ix0 = find(sx <= ux(u_idx)-alen,1,'last');
    if isempty(ix0) % Check edge cases (incomplete support for SAR aperture)
      ix0 = 1;
    end
    ix1 = find(sx >= ux(u_idx)+alen,1,'first');
    if isempty(ix1) % Check edge cases (incomplete support for SAR aperture)
      ix1 = length(sx);
    end

    anum = ix1-ix0+1;
    s_idx = zeros(1,ix1-ix0+1);
    tau = zeros(1,ix1-ix0+1);

    %loop over image z position
    large_search = true;
    % Compute the max-range pixel's WGS-84 elevation
%     uz_max = hz(u_idx) - (min(time(end),st(u_idx))*c/2 + max(0,time(end)-st(u_idx))*c/2/sqrt(er));
%     uz_max = hu(u_idx) - (min(time(end)*c/2,hu(u_idx)-su(u_idx)) + max(0,time(end)*c/2/sqrt(er)-(hu(u_idx)-su(u_idx))/sqrt(er)));
    for u_idz = start_range_bin:end_range_bin
      if time(u_idz) < 0
          continue;
        end
        % Print out status information
        if ~mod(u_idz-1,400)
          fprintf('  %d of %d, %d of %d (%s)\n', u_idx, nux, u_idz, nuz, datestr(now));
        end
        % Compute the current pixel's WGS-84 elevation
%         uz = hz(u_idx) - (min(time(u_idz),st(u_idx))*c/2 + max(0,time(u_idz)-st(u_idx))*c/2/sqrt(er));
          uz = hu(u_idx) - (min(time(u_idz)*c/2,hu(u_idx)-su(u_idx)) + max(0,time(u_idz)*c/2/sqrt(er)-(hu(u_idx)-su(u_idx))/sqrt(er)));
       
          theta_ice = atan((ux(u_idx)-sx(ix0:ix1))./(uz-sz(ix0:ix1)))+sm(ix0:ix1); %angle between the pixel and surface across the current aperture x-range.

        %loop over aperture (measurement location hz)
        data = zeros(ix1-ix0+1,1);
        for h_idx = ix0:1:ix1,
          if ~large_search && ((s_idx(h_idx-ix0+1) >= 1+11) && (s_idx(h_idx-ix0+1) <= anum-11)),
            % Otherwise use previous info and search locally (not sure how much this speeds it up)
            theta_air = atan((hx(h_idx)-sx(s_idx(h_idx-ix0+1)+ix0-1+(-10:10)))./(hz(h_idx)-sz(s_idx(h_idx-ix0+1)+ix0-1+(-10:10))))+sm(s_idx(h_idx-ix0+1)+ix0-1+(-10:10));
            [tmp,s_idx_temp] = min(abs(sin(theta_air)-sqrt(er)*sin(theta_ice(s_idx(h_idx-ix0+1)+(-10:10)))));
            s_idx(h_idx-ix0+1) = s_idx_temp + s_idx(h_idx-ix0+1) - 10 - 1;
          else
            % First time through and edge cases (-10 to 10): do a large search for snell equality
            large_search = false;
            theta_air = atan((hx(h_idx)-sx(ix0:ix1))./(hz(h_idx)-sz(ix0:ix1)))+sm(ix0:ix1); %angle between the antenna and surface across the current aperture x-range.
            [tmp,s_idx(h_idx-ix0+1)] = min(abs(sin(theta_air)-sqrt(er)*sin(theta_ice))); %minimum value satisfies snell
          end
          %calculate time delay
%           if ((uz+50) < su), %below the surface, use free space and ice
          if ((uz+50) < su(u_idx)), %below the surface, use free space and ice
            tau(h_idx-ix0+1) = (1/(c/2))*(sqrt((hx(h_idx)-sx(s_idx(h_idx-ix0+1)+ix0-1))^2 + (hz(h_idx)-sz(s_idx(h_idx-ix0+1)+ix0-1))^2) + ...
              sqrt(er)*sqrt((sx(s_idx(h_idx-ix0+1)+ix0-1)-ux(u_idx))^2 + (sz(s_idx(h_idx-ix0+1)+ix0-1)-uz)^2));
          else %above the surface, just use free space
            tau(h_idx-ix0+1) = (1/(c/2))*sqrt((hx(h_idx)-ux(u_idx))^2 + (hz(h_idx)-uz)^2);
          end
          %linear interpolation based on time delay
          iz = (tau(h_idx-ix0+1)-time(1))/dt + 1;
          if iz >= 1 && iz < nuz
            a1 = 1-mod(iz,1);
            a2 = 1-a1;
            data(h_idx-ix0+1) = a1*exp(1i*2*pi*fc*tau(h_idx-ix0+1)) * tx(floor(iz),h_idx) + ...
              a2*exp(1i*2*pi*fc*tau(h_idx-ix0+1)) * tx(floor(iz)+1,h_idx);
          end
        end
        
        % Remove phase associated with increasing depth (effect is to make
        %   SAR focussing compensate only for phase differences relative to
        %   the range of closest approach
        data = data.*exp(-1i*2*pi*fc*time(u_idz));
        
%         anump = round(((agl_eff+uz)/(agl_eff+uz_max))*anum); %scale aperture to depth *** need to include surface height in this calculation ***
        anump = round((500+(su(u_idx)-uz)/sqrt(er))/(500+1000/sqrt(er))*anum); % scale aperture to depth *** need to include surface height in this calculation ***
        if anump >anum
          anump = anum;
        end
        if anump < 20
          anump = 20;
        end
        blenp = 2*floor((anump-1)/(bins+1))+1; %calculte bin sizes
        anump = 0.5*(bins+1)*(blenp-1)+1; %tweak aperture size to accomodate all doppler bins
        a_off = floor(0.5*(anum-anump));
        awin = hanning(blenp);
        for idb = 1:bins, % coherently integrate all looks
            i0 = a_off + 1 + (idb-1)*(blenp-1)/2;
            i1 = i0 + blenp -1;
            u(u_idz,u_idx,idb) = mean(data(i0:i1).*awin);
        end;
    end;
end;
