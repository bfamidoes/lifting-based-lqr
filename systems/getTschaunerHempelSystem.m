function [A, B, orbit] = getTschaunerHempelSystem(p, a, e, f0, uRef)
% getTschaunerHempelSystem - Euler-discretized Tschauner-Hempel model, uniform steps in true anomaly
%
% State in scaled TH coordinates (see buildLvlh2Th); input the LVLH acceleration in units of uRef.
%
% Inputs:
%    p    - Steps per orbit
%    a, e - Chief semi-major axis (m) and eccentricity
%    f0   - True anomaly at step 0 (rad)
%    uRef - Input scaling (km/s^2)
%
% Outputs:
%    A, B  - State and input matrices, cell arrays
%    orbit - Orbit data (mu, h, n, e, f0, df, uRef), in km and s

mu = 3.986004418e5;                 % Earth gravitational parameter
aKm = a / 1e3;
h = sqrt(mu * aKm * (1 - e^2));     % specific angular momentum, km^2/s
df = 2 * pi / p;

A = cell(1, p);
B = cell(1, p);
for l = 1:p
    [Ac, Bc] = getThContinuousDynamics(f0 + (l - 1) * df, e, mu, h);
    A{l} = eye(6) + df * Ac;
    B{l} = df * Bc * uRef;
end

orbit = struct('mu', mu, 'h', h, 'n', sqrt(mu / aKm^3), 'e', e, 'f0', f0, 'df', df, 'uRef', uRef);

end
