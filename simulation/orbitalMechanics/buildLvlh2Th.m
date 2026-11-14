function T = buildLvlh2Th(f, e, mu, h)
% buildLvlh2Th - Transformation from LVLH to scaled Tschauner-Hempel coordinates
%
% Maps the LVLH state [position (km); velocity (km/s)] to TH coordinates, whose derivatives are
% taken with respect to the true anomaly.
%
% Inputs:
%    f  - True anomaly (rad)
%    e  - Eccentricity
%    mu - Gravitational parameter (km^3/s^2)
%    h  - Specific angular momentum (km^2/s)
%
% Outputs:
%    T - Transformation matrix (6-by-6)

rho = 1 + e*cos(f);

T = [mu*rho/h^(3/2)*eye(3),       zeros(3);
     -mu*e*sin(f)/h^(3/2)*eye(3),  h^(3/2)/(mu*rho)*eye(3)];

end
