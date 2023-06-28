function [theta,Ra,Ri,err]=refraction(h,d,yn,surfaceAngle,Eta_ice)
% [theta,Ra,Ri]=refraction(h,d,yn,surfaceAngle,Eta_ice)
% Function to calculate the refraction at the air-ice interface
% Linear interface assumed 
% h = airplane height in meters
% d = target depth in meters
% yn = alont-track offset (in meters)of the target from the atenna phase centers 
% surfaceAngle = ice surface slope angle in radians
% Eta_ice = air-ice refraction index
% theta = incidence angle in radians, the angle between the vertical and the propogation path in air
% Ra = distance from antenna phase center to air-ice interface
% Ri = distance from air-ice interface to target
% err: the mismatch of snell's law equation sin(theta)/sin(theta_i)-Eta_ice; 

% Jilu Li, Anthony, November 2012
% 

if Eta_ice == 1
    theta = atan(yn/(h+d));
    Ra = h/cos(theta);
    Ri = sqrt((h+d)^2+yn^2)-Ra;
    err = 0;
    return
end

% surface slope correction
% ==== modified to calculate rotation from bi-sector point
s           = sqrt( (h+d).^2 + yn.^2 );
newAngle    = 0.5*pi - asin(yn/s) + surfaceAngle;
h = h * cos(surfaceAngle);     % new h
sAbove      = h / sin(newAngle);
sBelow      = s - sAbove;
d           = sBelow * sin (newAngle);  % new d
if abs(newAngle-0.5*pi)<1e-6
    Ra = h;
    Ri = d;
    theta =0;
    err = 0;
else
yn          = sqrt( s.^2 - (h+d).^2 );  % new yn
% ====
    c4=Eta_ice^2-1;
    c3=-2*yn*c4;
    c2=c4*yn^2+(Eta_ice*h)^2-d^2;
    c1=-2*yn*(Eta_ice*h)^2;
    c0=(Eta_ice*yn*h)^2;
    poly=[c4 c3 c2 c1 c0];
    ya=roots(poly);
    ya_idx = find(abs(imag(ya))<=1e-5 & real(ya)>=0 & real(ya)<=yn+1e-5);
    if isempty(ya_idx)
        Ra = [];
        Ri =[];
        theta = [];
        err = [];
        return
    end
    ya=ya(ya_idx);
    if size(ya,1)>1 | ~isreal(ya(1))
        ya=real(ya(1));
    end
    Ra=sqrt(h^2+ya^2);
    Ri=sqrt(d^2+(yn-ya)^2);
    theta=sign(pi/2-newAngle)*asin(ya/Ra);
    sin_theta_r = (yn-ya)/Ri;
    err =sin(theta)/sin_theta_r-Eta_ice;
end