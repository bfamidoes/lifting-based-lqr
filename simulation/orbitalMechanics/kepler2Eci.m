function [r, v] = kepler2Eci(a, e, inc, argp, raan, f, mu)
% kepler2Eci - ECI position and velocity from Keplerian elements
%
% Computes the state in the perifocal frame and rotates it into ECI.
%
% Inputs:
%    a               - Semi-major axis (m)
%    e               - Eccentricity
%    inc, argp, raan - Inclination, argument of perigee and RAAN (rad)
%    f               - True anomaly (rad)
%    mu              - Gravitational parameter (m^3/s^2)
%
% Outputs:
%    r, v - Position and velocity in ECI (m, m/s)

semiLatusRectum = a * (1 - e^2);
rPerifocal = semiLatusRectum / (1 + e * cos(f)) * [cos(f); sin(f); 0];
vPerifocal = sqrt(mu / semiLatusRectum) * [-sin(f); e + cos(f); 0];

rotZ = @(x) [cos(x), -sin(x), 0; sin(x), cos(x), 0; 0, 0, 1];
rotX = @(x) [1, 0, 0; 0, cos(x), -sin(x); 0, sin(x), cos(x)];
perifocalToEci = rotZ(raan) * rotX(inc) * rotZ(argp);

r = perifocalToEci * rPerifocal;
v = perifocalToEci * vPerifocal;

end
