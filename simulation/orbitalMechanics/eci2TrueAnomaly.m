function f = eci2TrueAnomaly(r, v, mu)
% eci2TrueAnomaly - True anomaly from ECI position and velocity
%
% Angle from the eccentricity vector to r, about the angular momentum; requires e > 0.
%
% Inputs:
%    r, v - Position and velocity in ECI (m, m/s)
%    mu   - Gravitational parameter (m^3/s^2)
%
% Outputs:
%    f - True anomaly (rad), in (-pi, pi]

h = cross(r, v);
eVec = cross(v, h) / mu - r / norm(r);  % Eccentricity Vector
f = atan2(dot(h, cross(eVec, r)) / norm(h), dot(eVec, r));

end
