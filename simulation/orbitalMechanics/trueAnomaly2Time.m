function t = trueAnomaly2Time(f, e, n)
% trueAnomaly2Time - Time since perigee from true anomaly
%
% Converts via the eccentric and mean anomalies; f may span several revolutions.
%
% Inputs:
%    f - True anomaly (rad)
%    e - Eccentricity
%    n - Mean motion (rad/s)
%
% Outputs:
%    t - Time since perigee (s)

nRevs = floor(f / (2*pi));
fWrap = f - nRevs * 2 * pi;

E = 2 * atan2(sqrt(1 - e) * sin(fWrap/2), sqrt(1 + e) * cos(fWrap/2));
M = E - e * sin(E) + nRevs * 2 * pi;

t = M / n;

end
